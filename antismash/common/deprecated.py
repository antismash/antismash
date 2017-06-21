# License: GNU Affero General Public License v3 or later
# A copy of GNU AGPL v3 should have been included in this software package in LICENSE.txt.

"""
This file will be removed as soon as the new abstraction layer is completed.
"""

import Bio
import logging
import os
import re
import sys

from argparse import Namespace
from Bio.Seq import Seq
from Bio.SeqFeature import SeqFeature, FeatureLocation # for others importing
from Bio.SeqRecord import SeqRecord
from collections import defaultdict
from helperlibs.bio import seqio

from antismash.common import gff_parser


# temporary code skip logging # TODO
import inspect
import linecache

def CODE_SKIP_WARNING():
    prev = inspect.currentframe().f_back
    logging.critical("skipping code:" + prev.f_code.co_name +"():"+ linecache.getline(prev.f_code.co_filename, prev.f_lineno + 1))
# end temp


def get_cluster_features(seq_record):
    "Return all cluster features for a seq_record"
    return get_all_features_of_type(seq_record, "cluster")

def get_all_features_of_type(seq_record, types):
    "Return all features of the specified types for a seq_record"
    if isinstance(types, str):
        # force into a tuple
        types = (types, )
    features = []
    for f in seq_record.features:
        if f.type in types:
            features.append(f)
    return features

def get_cds_features(seq_record):
    "Return all CDS features for a seq_record"
    return get_all_features_of_type(seq_record, "CDS")

def get_withincluster_cds_features(seq_record):
    features = get_cds_features(seq_record)
    clusters = get_cluster_features(seq_record)
    withinclusterfeatures = []
    for feature in features:
        for cluster in clusters:
            if not (cluster.location.start <= feature.location.start <= cluster.location.end or \
               cluster.location.start <= feature.location.end <= cluster.location.end):
                continue
            if feature not in withinclusterfeatures:
                withinclusterfeatures.append(feature)
    return withinclusterfeatures
    
    
def parse_input_sequence(filename, options, genefinding):
    "Parse the input sequences from given filename"
    logging.info('Parsing input sequence %r', filename)

    sequences = []
    if not os.path.exists(filename):
        msg = "Sequence file not found: %r" % filename
        logging.error(msg)
        raise ValueError(msg)

    try:
        record_list = list(seqio.parse(filename))
        if len(record_list) == 0:
            logging.error('No sequence in file %r', filename)
        sequences.extend([rec for rec in record_list if len(rec.seq) > options.minlength or \
            ('contig' in rec.annotations or 'wgs_scafld' in rec.annotations or \
            'wgs' in rec.annotations)])
    except (ValueError, AssertionError) as e:
        logging.error('Parsing %r failed: %s', filename, e)
        sys.exit(1)
    except Exception as e:
        logging.error('Parsing %r failed with unhandled exception: %s',
                      filename, e)
        raise
    #Check if seq_records have appropriate content
    i = 0
    while i < len(sequences):
        sequence = sequences[i]
        sequence.seq = Seq(str(sequence.seq).replace("-", "").replace(":", ""))
        #Check if seq_record has either a sequence or has at least 80% of CDS features with 'translation' qualifier
        cdsfeatures = get_cds_features(sequence)
        cdsfeatures_with_translations = sum([1 for cdsfeature in cdsfeatures if 'translation' in cdsfeature.qualifiers])
        if len(sequence.seq) == 0 or (
                options.input_type == 'nucl' and \
                len(str(sequence.seq).replace("N","")) == 0 and \
                cdsfeatures_with_translations < 0.8 * len(cdsfeatures)):
            logging.error("Record %s has no sequence, skipping.", sequence.id)
            sequences.pop(i)
            continue

        if options.input_type == 'prot':
            if is_nucl_seq(sequence.seq):
                logging.error("Record %s is a nucleotide record, skipping.", sequence.id)
                sequences.pop(i)
                continue
        elif options.input_type == 'nucl':
            if not isinstance(sequence.seq.alphabet, Bio.Alphabet.NucleotideAlphabet) and not is_nucl_seq(sequence.seq):
                logging.error("Record %s is a protein record, skipping.", sequence.id)
                sequences.pop(i)
                continue
            sequence.seq.alphabet = Bio.Alphabet.generic_dna

        i += 1

    #If protein input, convert all protein seq_records to one nucleotide seq_record
    if options.input_type == 'prot':
        sequences = generate_nucl_seq_record(sequences)

    #Handle WGS master or supercontig entries
    check_for_wgs_scaffolds(sequences)

    #Now remove small contigs < minimum length again
    sequences = [rec for rec in sequences if len(rec.seq) > options.minlength]

    # Make sure we don't waste weeks of runtime on huge records, unless requested by the user
    old_len = len(sequences)
    if options.limit > -1:
        sequences = sequences[:options.limit]

    new_len = len(sequences)
    if new_len < old_len:
        options.triggered_limit = True
        logging.warning("Only analysing the first %d records (increase via --limit)", options.limit)

    #Check if no duplicate locus tags / gene IDs are found
    check_duplicate_gene_ids(sequences)

    #Store IDs for all entries
    CODE_SKIP_WARNING()
#    options.all_record_ids = {seq.id for seq in sequences}

    # For retaining the correct contig numbers, a second counter is required also including removed sequences without genes
    CODE_SKIP_WARNING()
#    options.orig_record_idx = 1

    # Check GFF suitability
    if options.genefinding_gff3:
        gff_parser.check_gff_suitability(options, sequences)

    i = 0
    while i < len(sequences):
        sequence = sequences[i]
        #Fix sequence name (will be ID) if it contains illegal chars
        illegal_chars  = '''!"#$%&()*+,:; \r\n\t=>?@[]^`'{|}/ '''
        for char in sequence.name:
            if char in illegal_chars:
                sequence.name = sequence.name.replace(char, "_")
        #Iterate through sequence objects
        if len(get_cds_features(sequence)) < 1:
            if options.gff3:
                logging.info("No CDS features found in record %r but GFF3 file provided, running GFF parser.", sequence.id)
                gff_parser.run(sequence, options)
                check_duplicate_gene_ids(sequences)
            else:
                logging.info("No CDS features found in record %r, running gene finding.", sequence.id)
                genefinding.run_on_record(sequence, options)
            if len(get_cds_features(sequence)) < 1:
                logging.info("No genes found, skipping record")
                sequences.pop(i)
                CODE_SKIP_WARNING()
#                options.orig_record_idx += 1
                continue
        #Fix locus tags
        fix_locus_tags(sequence)
        options.next_record_index()
        CODE_SKIP_WARNING()
#        options.orig_record_idx += 1
        i += 1

    #Make sure that all CDS entries in all seq_records have translation tags, otherwise add them
    add_translations(sequences)

    #Make sure that all seq_records have a sequence
    add_seq_record_seq(sequences)

    if len(sequences) > 1:
        options.start = -1
        options.end = -1
        logging.info("Discarding --start and --end options, as multiple entries are used.")

    i = 0
    while i < len(sequences):
        sequence = sequences[i]

        if options.start > 1:
            if options.start > len(sequence):
                logging.error('Specified analysis start point is at %r, which is larger ' \
                              'than record size %r', options.start, len(sequence))
                sys.exit(1)
            sequence = sequence[options.start-1:]
            # new sequence is shorter, so fix the end calculation
            options.end -= options.start
            sequences[i] = sequence

        if options.end > 0:
            if options.end > len(sequence):
                logging.error('Specified analysis end point is at %r, which is larger ' \
                              'than record size %r', options.end, len(sequence))
                sys.exit(1)
            sequence = sequence[:options.end]
            sequences[i] = sequence

        # Some programs write gaps as - not N, but translate() hates that
        if sequence.seq.find('-') > -1:
            sequence.seq = Seq(str(sequence.seq).replace('-', 'N'),
                               alphabet=sequence.seq.alphabet)

        # Some programs like to write gaps as X, translate() hates that
        if sequence.seq.find('X') > -1:
            sequence.seq = Seq(str(sequence.seq).replace('X', 'N'),
                               alphabet=sequence.seq.alphabet)
        if sequence.seq.find('x') > -1:
            sequence.seq = Seq(str(sequence.seq).replace('x', 'N'),
                               alphabet=sequence.seq.alphabet)

        i += 1

    #Fix sequence record IDs to be unique
    ids_used = []
    for sequence in sequences:
        seq_id = sequence.id
        if seq_id not in ids_used:
            ids_used.append(seq_id)
        else:
            x = 0
            #Make sure the length of the ID does not exceed 16
            if len(seq_id) <= 12:
                while "%s_%i" % (seq_id, x) in ids_used:
                    x += 1
                sequence.id = "%s_%i" % (seq_id, x)
            else:
                while "%s_%i" % (seq_id[:-4], x) in ids_used:
                    x += 1
                sequence.id = "%s_%i" % (seq_id[:-4], x)
            ids_used.append(sequence.id)
            options.all_record_ids.add(sequence.id) #Update all_record_ids with new record
    return sequences

def is_nucl_seq(sequence):
    other = str(sequence).lower()
    for char in "acgtn":
        other = other.replace(char, "")
    return len(other) < 0.2 * len(sequence)

def generate_nucl_seq_record(sequences):
    "Generate nucleotide seq_record"
    if len(sequences) == 0:
        return []
    seq_record = SeqRecord(Seq(""),id="Protein_Input", name="ProteinInput",
                   description="antiSMASH protein input")
    position = 0
    cds_features = []
    cdsnames = []
    for sequence in sequences:
        startpos = position
        endpos = position + len(sequence) * 3
        position += len(sequence) * 3 + 1000
        location = FeatureLocation(startpos, endpos)
        cdsfeature = SeqFeature(location, type="CDS")
        cdsfeature.strand = 1
        sequence_id = sequence.id[:15].replace(" ","_")
        if sequence_id not in cdsnames:
            cdsfeature.qualifiers['product'] = [sequence_id]
            cdsfeature.qualifiers['locus_tag'] = [sequence_id]
            cdsnames.append(sequence_id)
        else:
            x = 1
            while sequence_id[:8] + "_" + str(x) in cdsnames:
                x += 1
            cdsfeature.qualifiers['product'] = [sequence_id[:8] + "_" + str(x)]
            cdsfeature.qualifiers['locus_tag'] = [sequence_id[:8] + "_" + str(x)]
            cdsnames.append(sequence_id[:8] + "_" + str(x))
        cdsfeature.qualifiers['translation'] = [str(sequence.seq).replace('.', 'X')]
        cds_features.append(cdsfeature)
    seq_record.features.extend(cds_features)
    return [seq_record]

def check_duplicate_gene_ids(sequences):
    "Fix duplicate locus tags so that they are different"
    NO_TAG = "no_tag_found"
    high_water_mark = 0
    all_ids = defaultdict(lambda: False)
    for sequence in sequences:
        seq_ids = get_cds_features(sequence)
        for cdsfeature in seq_ids:
            gene_id = get_gene_id(cdsfeature)
            if not all_ids[gene_id]:
                all_ids[gene_id] = True
            else:
                if gene_id == NO_TAG:
                    x = high_water_mark + 1
                else:
                    x = 1
                id_str = "%s_%s" % ( gene_id[:8], x)
                while all_ids[id_str]:
                    x += 1
                    id_str = "%s_%s" % ( gene_id[:8], x)
                logging.debug("generated id %r", id_str)
                cdsfeature.qualifiers['product'] = [id_str]
                cdsfeature.qualifiers['locus_tag'] = [id_str]
                all_ids[id_str] = True
                if gene_id == NO_TAG:
                    high_water_mark = x


def fix_locus_tags(seq_record):
    "Fix CDS feature that don't have a locus_tag, gene name or protein id"
    next_locus_tag = 1

    cds_list = get_cds_features(seq_record)
    for feature in cds_list:
        if get_gene_id(feature) == "no_tag_found":
            feature.qualifiers['locus_tag'] = ['AUTOORF_%05d' % next_locus_tag]
            next_locus_tag += 1
        #Fix locus tags, gene names or protein IDs if they contain illegal chars
        illegal_chars  = '''!"#$%&()*+,:; \r\n\t=>?@[]^`'{|}/ '''
        if 'locus_tag' in feature.qualifiers:
            for char in feature.qualifiers['locus_tag'][0]:
                if char in illegal_chars:
                    feature.qualifiers['locus_tag'][0] = feature.qualifiers['locus_tag'][0].replace(char, "_")
        if 'gene' in feature.qualifiers:
            for char in feature.qualifiers['gene'][0]:
                if char in illegal_chars:
                    feature.qualifiers['gene'][0] = feature.qualifiers['gene'][0].replace(char, "_")
        if 'protein_id' in feature.qualifiers:
            for char in feature.qualifiers['protein_id'][0]:
                if char in illegal_chars:
                    feature.qualifiers['protein_id'][0] = feature.qualifiers['protein_id'][0].replace(char, "_")
                    
def add_translations(seq_records):
    "Add a translation qualifier to all CDS features"
    for seq_record in seq_records:
        logging.debug("Adding translations to record: %s", seq_record.id)
        cdsfeatures = get_cds_features(seq_record)
        for cdsfeature in cdsfeatures:
            if cdsfeature.qualifiers.get('translation'):
                continue
            if len(seq_record.seq) == 0:
                logging.error('No amino acid sequence in input entry for CDS %r, ' \
                        'and no nucleotide sequence provided to translate it from.', cdsfeature.id)
                raise ValueError("Missing sequence info for CDS %r" % cdsfeature.id)
            try:
                translation = str(get_aa_translation(seq_record, cdsfeature))
            except Bio.Data.CodonTable.TranslationError as e:
                logging.error('Getting amino acid sequences from %s, CDS %r failed: %s',
                        seq_record.name, cdsfeature.id, e)
                raise
            cdsfeature.qualifiers['translation'] = [translation]

def get_gene_id(feature):
    "Get the gene ID from locus_tag, gene name or protein id, in that order"
    if 'locus_tag' in feature.qualifiers:
        return feature.qualifiers['locus_tag'][0]
    if 'gene' in feature.qualifiers:
        return feature.qualifiers['gene'][0]
    if 'protein_id' in feature.qualifiers:
        return feature.qualifiers['protein_id'][0]
    return "no_tag_found"

def add_seq_record_seq(seq_records):
    for seq_record in seq_records:
        if len(seq_record.seq) == 0:
            seqmax = max([cds.location.start for cds in get_cds_features(seq_record)] + [cds.location.end for cds in get_cds_features(seq_record)])
            seq_record.seq = Seq(seqmax * "n")

def get_aa_translation(seq_record, feature):
    """Obtain content for translation qualifier for specific CDS feature in sequence record"""
    fasta_seq = feature.extract(seq_record.seq).ungap('-').translate(to_stop=True)
    if len(fasta_seq) == 0:
        logging.debug("Retranslating %s with stop codons", feature.id)
        fasta_seq = feature.extract(seq_record.seq).ungap('-').translate()
    if "*" in str(fasta_seq):
        fasta_seq = Seq(str(fasta_seq).replace("*","X"), Bio.Alphabet.generic_protein)
    if "-" in str(fasta_seq):
        fasta_seq = Seq(str(fasta_seq).replace("-",""), Bio.Alphabet.generic_protein)

    return fasta_seq


def check_for_wgs_scaffolds(seq_records):
    for seq_record in seq_records:
        #Check if seq_record is a WGS master record or a supercontig record
        if 'wgs_scafld' in seq_record.annotations \
                or 'wgs' in seq_record.annotations \
                or 'contig' in seq_record.annotations:
            raise RuntimeError("Incomplete whole genome shotgun records not supported")

def get_feature_dict(seq_record):
    """Get a dictionary mapping features to their IDs"""
    features = get_cds_features(seq_record)
    feature_by_id = {}
    for feature in features:
        gene_id = get_gene_id(feature)
        feature_by_id[gene_id] = feature
    return feature_by_id


def get_multifasta(seq_record):
    """Extract multi-protein FASTA from all CDS features in sequence record"""
    features = get_cds_features(seq_record)
    all_fastas = []
    for feature in features:
        gene_id = get_gene_id(feature)
        fasta_seq = feature.qualifiers.get('translation', [''])[0]
        if "-" in str(fasta_seq):
            fasta_seq = Seq(str(fasta_seq).replace("-",""), Bio.Alphabet.generic_protein)

        # Never write empty fasta entries
        if len(fasta_seq) == 0:
            logging.debug("No translation for %s, skipping", gene_id)
            assert feature.type.lower() != "cds"
            continue

        all_fastas.append(">%s\n%s" % (gene_id, fasta_seq))
    full_fasta = "\n".join(all_fastas)
    return full_fasta

def get_cluster_type(cluster):
    "Get product type of a gene cluster"
    return cluster.qualifiers['product'][0]

def get_cluster_cds_features(cluster, seq_record):
    clustercdsfeatures = []
    for feature in seq_record.features:
        if feature.type != 'CDS':
            continue
        if cluster.location.start <= feature.location.start <= cluster.location.end or \
           cluster.location.start <= feature.location.end <= cluster.location.end:
            clustercdsfeatures.append(feature)
    return clustercdsfeatures

def strip_record(seq_record):

    new_features = []

    for feature in seq_record.features:

        # Discard features added by antiSMASH
        if feature.type in ('cluster', 'cluster_border', 'CDS_motif', 'aSDomain'):
            continue

        # clean up antiSMASH annotations in CDS features
        if feature.type == 'CDS':
            if 'sec_met' in feature.qualifiers:
                del feature.qualifiers['sec_met']

        new_features.append(feature)

    seq_record.features = new_features


def sort_features(seq_record):
    "Sort features in a seq_record by their position"

    #Check if all features have a proper location assigned
    for feature in seq_record.features:
        if feature.location is None:
            if feature.id != "<unknown id>":
                logging.error("Feature '%s' has no proper location assigned", feature.id)
            elif "locus_tag" in feature.qualifiers:
                logging.error("Feature '%s' has no proper location assigned", feature.qualifiers["locus_tag"][0])
            else:
                logging.error("File contains feature without proper location assignment")
            raise ValueError("File contains feature without proper location assignment")
    #Sort features by location
    seq_record.features.sort(key=lambda x: (x.location.start, x.location.end))

def fix_record_name_id(seq_record, options):
    "Fix a seq record's name and id to be <= 16 characters, the GenBank limit; if record name is too long, add c000X prefix"

    def _shorten_ids(idstring, options):
        contigstrmatch = re.search(r"onti?g?(\d+)\b", idstring)
        if not contigstrmatch:
            # if there is a substring "[Ss]caf(fold)XXX" use this number
            contigstrmatch = re.search(r"caff?o?l?d?(\d+)\b", idstring)
        if not contigstrmatch:
            # if there is a substring " cXXX" use this number
            contigstrmatch = re.search(r"\bc(\d+)\b", idstring)
        if contigstrmatch:
            contig_no = int(contigstrmatch.group(1))
        else:
            # if the contig number cannot be parsed out, just count the contigs from 1 to n
            contig_no = options.orig_record_idx

        return "c{ctg:05d}_{origid}..".format(ctg=contig_no, origid=idstring[:7])

    if seq_record.id == "unknown.1":
        seq_record.id = "unk_seq_{ctg:05d}".format(ctg=options.orig_record_idx)
        logging.warn('Invalid sequence id "unknown.1", replaced by %s', seq_record.id)

    if seq_record.name == "unknown":
        seq_record.name = "unk_seq_{ctg:05d}".format(ctg=options.orig_record_idx)
        logging.warn('Invalid sequence name "unknown", replaced by %s', seq_record.name)

    if len(seq_record.id) > 16:
        oldid = seq_record.id

        #Check if it is a RefSeq accession number like NZ_AMZN01000079.1 that is just too long because of the version number behind the dot
        if (seq_record.id[-2] == "." and
                seq_record.id.count(".") == 1 and
                len(seq_record.id.partition(".")[0]) <= 16 and
                seq_record.id.partition(".")[0] not in options.all_record_ids):
            seq_record.id = seq_record.id.partition(".")[0]
            options.all_record_ids.add(seq_record.id)
        else: #Check if the ID suggested by _shorten_ids is unique
            if _shorten_ids(oldid, options) not in options.all_record_ids:
                seq_record.id = _shorten_ids(oldid, options)
                options.all_record_ids.add(seq_record.id)
            else:
                x = 0
                while "%s_%i" % (seq_record.id[:16][:-4], x) in options.all_record_ids:
                    x += 1
                seq_record.id = "%s_%i" % (seq_record.id[:16][:-4], x)
                options.all_record_ids.add(seq_record.id)

        logging.warn('Fasta header too long: renamed "%s" to "%s"', oldid, seq_record.id)
        if seq_record.id not in options.extrarecord:
            options.extrarecord[seq_record.id] = Namespace()
        if "extradata" not in options.extrarecord[seq_record.id]:
            options.extrarecord[seq_record.id].extradata = {}
        if "orig_id" not in options.extrarecord[seq_record.id].extradata:
            options.extrarecord[seq_record.id].extradata["orig_id"] = oldid

    if len(seq_record.name) > 16:

        seq_record.name = _shorten_ids(seq_record.name, options)

    if 'accession' in seq_record.annotations and \
       len(seq_record.annotations['accession']) > 16:
        acc = seq_record.annotations['accession']

        seq_record.annotations['accession'] = _shorten_ids(acc, options)

    # Remove illegal characters from name: otherwise, file cannot be written
    illegal_chars  = '''!"#$%&()*+,:;=>?@[]^`'{|}/ '''
    for char in seq_record.id:
        if char in illegal_chars:
            seq_record.id = seq_record.id.replace(char,"")
    for char in seq_record.name:
        if char in illegal_chars:
            seq_record.id = seq_record.name.replace(char,"")

