# License: GNU Affero General Public License v3 or later
# A copy of GNU AGPL v3 should have been included in this software package in LICENSE.txt.

"""
This file will be removed as soon as all modules from antiSMASH 4 have been
converted
"""

import logging
import os
import re
import sys

import Bio
from Bio.Data.IUPACData import protein_letters
from Bio.Seq import Seq
from Bio.SeqFeature import SeqFeature, FeatureLocation # for others importing
from Bio.SeqRecord import SeqRecord
from Bio.SeqUtils.ProtParam import ProteinAnalysis
from helperlibs.bio import seqio

from antismash.common import gff_parser
from antismash.common.all_orfs import scan_orfs, sort_orfs
from antismash.common.secmet import Record, CDSFeature, Feature

# temporary code skip logging # TODO
import inspect
import linecache

def CODE_SKIP_WARNING():
    prev = inspect.currentframe().f_back
    logging.critical("skipping code:" + prev.f_code.co_name +"():" \
            + linecache.getline(prev.f_code.co_filename, prev.f_lineno + 1).replace('%', '%%'))
# end temp


def parse_input_sequence(filename, options):
    "Parse the input sequences from given filename"
    logging.info('Parsing input sequence %r', filename)

    sequences = []
    if not os.path.exists(filename):
        msg = "Sequence file not found: %r" % filename
        logging.error(msg)
        raise ValueError(msg)

    try:
        record_list = list(seqio.parse(filename))
        if not record_list:
            logging.error('No sequence in file %r', filename)
        sequences.extend([rec for rec in record_list if len(rec.seq) > options.minlength or \
            ('contig' in rec.annotations or 'wgs_scafld' in rec.annotations or \
            'wgs' in rec.annotations)])
    except (ValueError, AssertionError) as err:
        logging.error('Parsing %r failed: %s', filename, err)
        raise
    except Exception as err:
        logging.error('Parsing %r failed with unhandled exception: %s',
                      filename, err)
        raise
    return [Record.from_biopython(sequence) for sequence in sequences]

def pre_process_sequences(sequences, options, genefinding):
    # keep count of how many records matched filter
    matching_filter = 0

    # Check if seq_records have appropriate content
    for i, sequence in enumerate(sequences):
        if options.limit_to_record and options.limit_to_record != sequence.id:
            sequence.skip = "did not match filter: %s" % options.limit_to_record
        else:
            matching_filter += 1

        sequence.record_index = i
        sequence.seq = Seq(str(sequence.seq).replace("-", "").replace(":", ""))
        # Check if seq_record has either a sequence or has at least 80% of CDS features with 'translation' qualifier
        cdsfeatures = sequence.get_cds_features()
        cdsfeatures_with_translations = sum([1 for cdsfeature in cdsfeatures if cdsfeature.translation])
        if not sequence.seq or (
                options.input_type == 'nucl' and \
                not str(sequence.seq).replace("N", "") and \
                cdsfeatures_with_translations < 0.8 * len(cdsfeatures)):
            logging.error("Record %s has no sequence, skipping.", sequence.id)
            sequence.skip = "contains no sequence"
            continue

        if options.input_type == 'prot':
            if is_nucl_seq(sequence.seq):
                logging.error("Record %s is a nucleotide record, skipping.", sequence.id)
                sequence.skip = "nucleotide record in protein mode"
                continue
        elif options.input_type == 'nucl':
            if not isinstance(sequence.seq.alphabet, Bio.Alphabet.NucleotideAlphabet) and not is_nucl_seq(sequence.seq):
                logging.error("Record %s is a protein record, skipping.", sequence.id)
                sequence.skip = "protein record in nucleotide mode"
                continue
            sequence.seq.alphabet = Bio.Alphabet.generic_dna

    if options.limit_to_record:
        limit = options.limit_to_record
        if matching_filter == 0:
            logging.error("No sequences matched filter: %s", limit)
            raise ValueError("No sequences matched filter: %s" % limit)
        elif matching_filter != len(sequences):
            logging.info("Skipped %d sequences not matching filter: %s",
                         len(sequences) - matching_filter, limit)

    #If protein input, convert all protein seq_records to one nucleotide seq_record
    if options.input_type == 'prot':
        sequences = generate_nucl_seq_record(sequences)

    #Handle WGS master or supercontig entries
    check_for_wgs_scaffolds(sequences)

    #Now remove small contigs < minimum length again
    for sequence in sequences:
        if len(sequence.seq) < options.minlength:
            sequence.skip = "smaller than minimum length (%d)" % options.minlength

    # Make sure we don't waste weeks of runtime on huge records, unless requested by the user
    warned = False
    if options.limit > -1:
        meaningful = 0
        for sequence in sequences:
            if sequence.skip:
                continue
            meaningful += 1
            if meaningful > options.limit:
                if not warned:
                    logging.warning("Only analysing the first %d records (increase via --limit)", options.limit)
                sequence.skip = "skipping all but first {0} meaningful records (--limit {0}) ".format(options.limit)

    if warned:
        options.triggered_limit = True

    #Check if no duplicate locus tags / gene IDs are found
    check_duplicate_gene_ids(sequences)

    # Check GFF suitability
    if options.genefinding_gff3:
        single_entry = gff_parser.check_gff_suitability(options, sequences)

    # Ensure all records have valid names
    for seq_record in sequences:
        fix_record_name_id(seq_record, {seq.id for seq in sequences}, options)

    for sequence in sequences:
        if sequence.skip:
            continue
        #Fix sequence name (will be ID) if it contains illegal chars
        illegal_chars = '''!"#$%&()*+,:; \r\n\t=>?@[]^`'{|}/ '''
        for char in sequence.name:
            if char in illegal_chars:
                sequence.name = sequence.name.replace(char, "_")
        #Iterate through sequence objects
        if len(sequence.get_cds_features()) < 1:
            if options.genefinding_gff3:
                logging.info("No CDS features found in record %r but GFF3 file provided, running GFF parser.", sequence.id)
                gff_parser.run(sequence, single_entry, options)
                check_duplicate_gene_ids(sequences)
            elif options.genefinding_tool != "none":
                logging.info("No CDS features found in record %r, running gene finding.", sequence.id)
                genefinding.run_on_record(sequence, options)
            if len(sequence.get_cds_features()) < 1:
                logging.info("No genes found, skipping record")
                sequence.skip = "No genes found"
                continue
        #Fix locus tags
        fix_locus_tags(sequence)

    #Make sure that all CDS entries in all seq_records have translation tags, otherwise add them
    add_translations(sequences)

    #Make sure that all seq_records have a sequence
    add_seq_record_seq(sequences)

    if len(sequences) > 1 and (options.start != -1 or options.end != -1):
        options.start = -1
        options.end = -1
        logging.info("Discarding --start and --end options, as multiple entries are used.")

    for i, sequence in enumerate(sequences):
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

        # Some programs write gaps as - or X or x, translate requires N
        for char in ['-', 'X', 'x']:
            sequence.seq = Seq(str(sequence.seq).replace(char, 'N'),
                               alphabet=sequence.seq.alphabet)



    #Fix sequence record IDs to be unique
    ids_used = []
    for sequence in sequences:
        seq_id = sequence.id
        if seq_id not in ids_used:
            ids_used.append(seq_id)
            continue
        prefix = seq_id
        if len(prefix) > 11:
            prefix = prefix[:11]


        suffix = 0
        #Make sure the length of the ID does not exceed 16
        if len(seq_id) <= 11:
            while "%s_%i" % (seq_id, suffix) in ids_used:
                suffix += 1
            sequence.id = "%s_%i" % (seq_id, suffix)
        else:
            while "%s_%i" % (seq_id[:-4], suffix) in ids_used:
                suffix += 1
            sequence.id = "%s_%i" % (seq_id[:-4], suffix)
        ids_used.append(sequence.id)
    return sequences

def generate_unique_id(prefix, existing_ids, start=0, max_length=-1):
    """ Generate a identifier of the form prefix_num, e.g. seq_15.

        Args:
            prefix: The text portion of the name.
            existing_ids: The current identifiers to avoid collision with.
            start: An integer to start counting at (default: 0)
            max_length: The maximum length allowed for the identifier,
                        values less than 1 are considerd to be no limit.

        Returns:
            A tuple of the identifier generated and the value of the counter
                at the time the identifier was generated, e.g. ("seq_15", 15)

    """
    counter = int(start)
    existing_ids = set(existing_ids)

    format_string = "{}_{}".format(prefix, counter)
    name = format_string % counter
    while name in existing_ids:
        counter += 1
        name = format_string % counter
    if max_length > 0 and len(name) > max_length:
        raise RuntimeError("Could not generate unique id for %s after %d iterations" % (prefix, counter - start))
    return name, counter

def is_nucl_seq(sequence):
    other = str(sequence).lower()
    for char in "acgtn":
        other = other.replace(char, "")
    return len(other) < 0.2 * len(sequence)

def generate_nucl_seq_record(sequences):
    "Generate single nucleotide seq_record from supplied sequences"
    if not sequences:
        raise ValueError("Cannot generate nucleotide records of empty input")
    record = Record(Seq(""), id="Protein_Input", name="ProteinInput",
                   description="antiSMASH protein input")
    position = 0
    cdsnames = set()
    for sequence in sequences:
        startpos = position
        endpos = position + len(sequence) * 3
        position += len(sequence) * 3 + 1000
        location = FeatureLocation(startpos, endpos, strand=1)
        name = sequence.id[:15].replace(" ", "_")
        if name in cdsnames:
            name, _ = generate_unique_id(name[:8], cdsnames)
        cdsnames.add(name)
        translation = str(sequence.seq).replace('.', 'X')
        cdsfeature = CDSFeature(location, translation, product=name, locus_tag=name)
        record.add_cds_feature(cdsfeature)
    return record

def check_duplicate_gene_ids(sequences):
    "Fix duplicate locus tags so that they are different"
    no_tag = "no_tag_found"
    high_water_mark = 0
    all_ids = set()
    for sequence in sequences:
        for cdsfeature in sequence.get_cds_features():
            name = cdsfeature.get_name()
            if not name:
                name = no_tag
            if name == no_tag:
                name, high_water_mark = generate_unique_id(name[:8], all_ids,
                                                start=high_water_mark + 1)
            elif name in all_ids:
                name, _ = generate_unique_id(name[:8], all_ids, start=1)
            if cdsfeature.product is None:
                cdsfeature.product = name
            cdsfeature.locus_tag = name
            all_ids.add(name)


def fix_locus_tags(seq_record):
    "Fix CDS feature that don't have a locus_tag, gene name or protein id"
    next_locus_tag = 1

    for feature in seq_record.get_cds_features():
        if feature.get_name() == "no_tag_found":
            feature.locus_tag = '%s_%05d' % (seq_record.id, next_locus_tag)
            next_locus_tag += 1
        #Fix locus tags, gene names or protein IDs if they contain illegal chars
        illegal_chars = '''!"#$%&()*+,:; \r\n\t=>?@[]^`'{|}/ '''
        for attr in ["locus_tag", "gene", "protein_id"]:
            val = getattr(feature, attr)
            if not val:
                continue
            for char in val:
                if char in illegal_chars:
                    val = val.replace(char, "_")
            setattr(feature, attr, val)

def add_translations(seq_records):
    "Add a translation qualifier to all CDS features"
    for seq_record in seq_records:
        if seq_record.skip:
            continue
        logging.debug("Adding translations to record: %s", seq_record.id)
        cdsfeatures = seq_record.get_cds_features()
        for cdsfeature in cdsfeatures:
            if cdsfeature.translation:
                continue
            if not seq_record.seq:
                logging.error('No amino acid sequence in input entry for CDS %r, ' \
                        'and no nucleotide sequence provided to translate it from.', cdsfeature.id)
                raise ValueError("Missing sequence info for CDS %r" % cdsfeature.id)
            try:
                translation = str(seq_record.get_aa_translation_of_feature(cdsfeature))
            except Bio.Data.CodonTable.TranslationError as err:
                logging.error('Getting amino acid sequences from %s, CDS %r failed: %s',
                        seq_record.name, cdsfeature.id, err)
                raise
            cdsfeature.translation = translation

def add_seq_record_seq(seq_records):
    for seq_record in seq_records:
        if not seq_record.seq:
            cds_features = seq_record.get_cds_features()
            start_max = max([cds.location.start for cds in cds_features])
            end_max = max([cds.location.end for cds in cds_features])
            seq_record.seq = Seq(max([start_max, end_max]) * "n")

def check_for_wgs_scaffolds(seq_records):
    for seq_record in seq_records:
        #Check if seq_record is a WGS master record or a supercontig record
        if 'wgs_scafld' in seq_record.annotations \
                or 'wgs' in seq_record.annotations \
                or 'contig' in seq_record.annotations:
            raise RuntimeError("Incomplete whole genome shotgun records not supported")

def get_feature_dict(seq_record):
    """Get a dictionary mapping features to their IDs"""
    features = seq_record.get_cds_features()
    feature_by_id = {}
    for feature in features:
        gene_id = feature.get_name()
        feature_by_id[gene_id] = feature
    return feature_by_id


def get_multifasta(seq_record):
    """Extract multi-protein FASTA from all CDS features in sequence record"""
    features = seq_record.get_cds_features()
    all_fastas = []
    for feature in features:
        gene_id = feature.get_name()
        fasta_seq = feature.translation
        if "-" in str(fasta_seq):
            fasta_seq = Seq(str(fasta_seq).replace("-", ""), Bio.Alphabet.generic_protein)

        # Never write empty fasta entries
        if not fasta_seq:
            logging.error("No translation for CDS %s", gene_id)
            raise ValueError("No translation for CDS %s" % gene_id)

        all_fastas.append(">%s\n%s" % (gene_id, fasta_seq))
    full_fasta = "\n".join(all_fastas)
    return full_fasta

def writefasta(names, seqs, filename):
    "Write sequence to a file"
    e = 0
    f = len(names) - 1
    out_file = open(filename,"w")
    while e <= f:
        out_file.write(">")
        out_file.write(names[e])
        out_file.write("\n")
        out_file.write(seqs[e])
        out_file.write("\n")
        e += 1
    out_file.close()

def strip_record(seq_record):
    """ Discard antismash specific features and feature qualifiers """
    seq_record.clear_clusters()
    seq_record.clear_cluster_borders()
    seq_record.clear_cds_motifs()
    seq_record.clear_antismash_domains()

    # clean up antiSMASH annotations in CDS features
    for feature in seq_record.get_cds_features():
        feature.sec_met = None

def fix_record_name_id(seq_record, all_record_ids, options):
    "Fix a seq record's name and id to be <= 16 characters, the GenBank limit; if record name is too long, add c000X prefix"

    def _shorten_ids(idstring):
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
            contig_no = seq_record.record_index

        return "c{ctg:05d}_{origid}..".format(ctg=contig_no, origid=idstring[:7])

    if seq_record.id == "unknown.1":
        seq_record.id = "unk_seq_{ctg:05d}".format(ctg=seq_record.record_index)
        logging.warning('Invalid sequence id "unknown.1", replaced by %s', seq_record.id)

    if seq_record.name == "unknown":
        seq_record.name = "unk_seq_{ctg:05d}".format(ctg=options.record_index)
        logging.warning('Invalid sequence name "unknown", replaced by %s', seq_record.name)

    if len(seq_record.id) > 16:
        oldid = seq_record.id

        #Check if it is a RefSeq accession number like NZ_AMZN01000079.1 that is just too long because of the version number behind the dot
        if (seq_record.id[-2] == "." and
                seq_record.id.count(".") == 1 and
                len(seq_record.id.partition(".")[0]) <= 16 and
                seq_record.id.partition(".")[0] not in all_record_ids):
            seq_record.id = seq_record.id.partition(".")[0]
            all_record_ids.add(seq_record.id)
        else: #Check if the ID suggested by _shorten_ids is unique
            if _shorten_ids(oldid) not in all_record_ids:
                name = _shorten_ids(oldid)
            else:
                name, _ = generate_unique_id(seq_record.id[:12], all_record_ids, max_length=16)
            seq_record.id = name
            all_record_ids.add(name)

        logging.warning('Fasta header too long: renamed "%s" to "%s"', oldid, seq_record.id)

    if len(seq_record.name) > 16:

        seq_record.name = _shorten_ids(seq_record.name)

    if 'accession' in seq_record.annotations and \
       len(seq_record.annotations['accession']) > 16:
        acc = seq_record.annotations['accession']

        seq_record.annotations['accession'] = _shorten_ids(acc)

    # Remove illegal characters from name: otherwise, file cannot be written
    illegal_chars = '''!"#$%&()*+,:;=>?@[]^`'{|}/ '''
    for char in seq_record.id:
        if char in illegal_chars:
            seq_record.id = seq_record.id.replace(char, "")
    for char in seq_record.name:
        if char in illegal_chars:
            seq_record.name = seq_record.name.replace(char, "")

def get_feature_dict_protein_id(record):
    logging.critical("get_feature_dict_protein_id(record) called, did you mean record.get_cds_mapping()?")
    return record.get_cds_mapping()


def get_smcog_annotations(seq_record):
    logging.critical("get_smcog_annotations(): should use secmet for smCOG note")
    smcogdict = {}
    smcogdescriptions = {}
    for feature in seq_record.get_cds_features():
        for note in feature.notes:
            if "smCOG: " in note:
                smcogid = note.partition("smCOG: ")[2].partition(":")[0]
                smcog_descr = note.partition("smCOG: ")[2].partition(":")[2].partition("(Score:")[0]
                smcogdict[feature.get_name()] = smcogid
                smcogdescriptions[smcogid] = smcog_descr
    return smcogdict, smcogdescriptions

def get_pksnrps_cds_features(seq_record):
    logging.critical("skipping get_pksnrps_cds_features(), missing PKSNRPS module")
    return []
#    features = get_cds_features(seq_record)
#    pksnrpscoregenes = []
#    for feature in features:
#        if 'sec_met' in feature.qualifiers:
#            for annotation in feature.qualifiers['sec_met']:
#                if annotation.startswith('NRPS/PKS Domain:'):
#                    pksnrpscoregenes.append(feature)
#                    break
#    return pksnrpscoregenes

def get_nrpspks_domain_dict(seq_record):
    logging.critical("skipping get_pksnrps_domain_dict(), missing PKSNRPS module")
    return {}
#    domaindict = {}
#    features = get_cds_features(seq_record)
#    for feature in features:
#        domainlist = []
#        if 'sec_met' in feature.qualifiers:
#            domains = [qualifier for qualifier in feature.qualifiers['sec_met'] if "NRPS/PKS Domain: " in qualifier]
#            for domain in domains:
#                hit_id =  domain.partition("NRPS/PKS Domain: ")[2].partition(" (")[0]
#                domstart = domain.partition("NRPS/PKS Domain: ")[2].partition(" (")[2].partition("-")[0]
#                domend = domain.partition("NRPS/PKS Domain: ")[2].partition("). ")[0].rpartition("-")[2]
#                evalue = domain.partition("E-value: ")[2].partition(". Score:")[0]
#                bitscore = domain.partition("Score: ")[2].partition(";")[0]
#                domainlist.append([hit_id, int(domstart), int(domend), evalue, float(bitscore)])
#            if len(domainlist) > 0:
#                domaindict[feature.get_name()] = domainlist
#    return domaindict

def get_nrpspks_substr_spec_preds(seq_record):
    logging.critical("skipping get_nrpspks_substr_spec_preds(), missing PKSNRPS module")
    class Dummy:
        pass
    substr_spec_preds = Dummy()
    substr_spec_preds.consensuspreds = {}
    substr_spec_preds.nrps_svm_preds = {}
    substr_spec_preds.nrps_code_preds = {}
    substr_spec_preds.minowa_nrps_preds = {}
    substr_spec_preds.pks_code_preds = {}
    substr_spec_preds.minowa_pks_preds = {}
    substr_spec_preds.minowa_cal_preds = {}
    substr_spec_preds.kr_activity_preds = {}
    substr_spec_preds.kr_stereo_preds = {}
    return substr_spec_preds
#    features = get_cds_features(seq_record)
#    for feature in features:
#        nrat, nra, nrcal, nrkr = 0, 0, 0, 0
#        if 'sec_met' in feature.qualifiers:
#            domains = [qualifier for qualifier in feature.qualifiers['sec_met'] if "NRPS/PKS Domain: " in qualifier]
#            for domain in domains:
#                if "AMP-binding" in domain or "A-OX" in domain:
#                    nra += 1
#                    domainname = feature.get_name() + "_A" + str(nra)
#                    predictionstext = domain.partition("Substrate specificity predictions:")[2]
#                    nrps_svm_pred = predictionstext.partition(" (NRPSPredictor2 SVM)")[0]
#                    nrps_code_pred = predictionstext.partition(" (NRPSPredictor2 SVM), ")[2].partition(" (Stachelhaus code)")[0]
#                    minowa_nrps_pred = predictionstext.partition("(Stachelhaus code), ")[2].partition(" (Minowa)")[0]
#                    consensuspred = predictionstext.partition("(Minowa), ")[2].partition(" (consensus)")[0]
#                    substr_spec_preds.nrps_svm_preds[domainname] = nrps_svm_pred
#                    substr_spec_preds.nrps_code_preds[domainname] = nrps_code_pred
#                    substr_spec_preds.minowa_nrps_preds[domainname] = minowa_nrps_pred
#                    substr_spec_preds.consensuspreds[domainname] = consensuspred
#                elif "PKS_AT" in domain:
#                    nrat += 1
#                    domainname = feature.get_name() + "_AT" + str(nrat)
#                    predictionstext = domain.partition("Substrate specificity predictions:")[2]
#                    pks_code_pred = predictionstext.partition(" (PKS signature)")[0]
#                    minowa_pks_pred = predictionstext.partition("(PKS signature), ")[2].partition(" (Minowa)")[0]
#                    consensuspred = predictionstext.partition("(Minowa), ")[2].partition(" (consensus)")[0]
#                    substr_spec_preds.pks_code_preds[domainname] = pks_code_pred
#                    substr_spec_preds.minowa_pks_preds[domainname] = minowa_pks_pred
#                    substr_spec_preds.consensuspreds[domainname] = consensuspred
#                elif "CAL_domain" in domain:
#                    nrcal += 1
#                    domainname = feature.get_name() + "_CAL" + str(nrcal)
#                    predictionstext = domain.partition("Substrate specificity predictions:")[2]
#                    minowa_cal_pred = predictionstext.partition(" (Minowa)")[0]
#                    consensuspred = predictionstext.partition("(Minowa), ")[2].partition(" (consensus)")[0]
#                    substr_spec_preds.minowa_cal_preds[domainname] = minowa_cal_pred
#                    substr_spec_preds.consensuspreds[domainname] = consensuspred
#                elif "PKS_KR" in domain:
#                    nrkr += 1
#                    domainname = feature.get_name() + "_KR" + str(nrkr)
#                    activityprediction = domain.partition("Predicted KR activity: ")[2].partition(";")[0]
#                    stereoprediction = domain.partition("Predicted KR stereochemistry: ")[2].partition(";")[0]
#                    substr_spec_preds.kr_activity_preds[domainname] = activityprediction
#                    substr_spec_preds.kr_stereo_preds[domainname] = stereoprediction
#    return substr_spec_preds

def get_structure_pred(cluster):
    "Return all a structure prediction for a cluster feature"
    for note in cluster.notes:
        if "Monomers prediction: " in note:
            return note.partition("Monomers prediction: ")[2]
    if cluster.get_product_string() == 'ectoine': # TODO: wat
        return 'ectoine'
    return "N/A"

def get_version():
    logging.critical("dummy get_version() being called")
    return "antismash-5.alpha"

class RobustProteinAnalysis(ProteinAnalysis):
    PROTEIN_LETTERS = set(protein_letters)
    def __init__(self, prot_sequence, monoisotopic=False, invalid="ignore"):
        if invalid not in ("ignore", "average"):
            raise ValueError("Invalid needs to be 'ignore' or 'average', not {!r}".format(invalid))
        self._invalid = invalid

        prot_sequence = prot_sequence.upper()

        self.original_sequence = prot_sequence
        # replace all invalids with ' '
        prot_sequence = []
        for i in self.original_sequence:
            if i in RobustProteinAnalysis.PROTEIN_LETTERS:
                prot_sequence.append(i)
        prot_sequence = "".join(prot_sequence)
        super(RobustProteinAnalysis, self).__init__(prot_sequence, monoisotopic)

    def molecular_weight(self):
        mw = super(RobustProteinAnalysis, self).molecular_weight()
        if self._invalid == "average":
            aa_difference = len(self.original_sequence) - len(self.sequence)
            mw += 110 * aa_difference

        return mw

def find_all_orfs(seq_record, cluster): # the old lassopeptides.find_all_orfs
    """Find all ORFs in gene cluster outside annotated CDS features"""
    # Get sequence just for the gene cluster
    fasta_seq = seq_record.seq[cluster.location.start:cluster.location.end]

    # Find orfs throughout the cluster
    forward_matches = scan_orfs(fasta_seq, 1, cluster.location.start)
    reverse_matches = scan_orfs(fasta_seq.complement(), -1, cluster.location.start)
    all_orfs = forward_matches + reverse_matches

    orfnr = 1
    new_features = []

    for orf in sort_orfs(all_orfs):
        # Remove if overlaps with existing CDSs
        skip = False
        for cds in cluster.cds_children:
            if orf.start in cds.location or orf.end in cds.location or cds.location.start in orf or cds.location.end in orf:
                skip = True
                break
        if skip:
            continue
        loc = orf
        dummy_feature = Feature(loc, feature_type="dummy")
        feature = CDSFeature(loc, str(seq_record.get_aa_translation_of_feature(dummy_feature)),
                             locus_tag='cluster_%s_allorf%03d' % (cluster.get_cluster_number(), orfnr))
        new_features.append(feature)
        orfnr += 1

    return new_features

def distance_to_pfam(seq_record, query, hmmer_profiles): #also from lassopeptides
    """Function to check how many nt a gene is away from a gene with one of a list of given Pfams"""
    nt = 40000 #maximum number of nucleotides distance to search
    #Get all CDS features in seq_record
    cds_features = seq_record.get_cds_features()
    #Get all CDS features within <X nt distances
    close_cds_features = []
    distance = {}
    for cds in cds_features:
        if query.location.start - nt <= cds.location.start <= query.location.end + nt or \
           query.location.start - nt <= cds.location.end <= query.location.end + nt:
            close_cds_features.append(cds)
            distance[cds.get_name()] = min([
                                abs(cds.location.start - query.location.end),
                                abs(cds.location.end - query.location.start),
                                abs(cds.location.start - query.location.start),
                                abs(cds.location.end - query.location.end)])
    #For nearby CDS features, check if they have hits to the pHMM
    closest_distance = -1
    for cds in close_cds_features:
        if cds.sec_met:
            for profile in hmmer_profiles:
                if profile in cds.sec_met.domains:
                    if closest_distance == -1 or distance[cds.get_name()] < closest_distance:
                        closest_distance = distance[cds.get_name()]
    return closest_distance

def get_specific_multifasta(features):
    """Extract multi-protein FASTA from provided features"""
    all_fastas = []
    for feature in features:
        all_fastas.append(">%s\n%s" % (feature.get_name(), feature.translation))
    return "\n".join(all_fastas)

def hmmlengths(hmmfile):
    lengths = {}
    with open(hmmfile,"r") as handle:
        contents = handle.read()
    contents = contents.replace("\r", "\n")
    hmms = contents.split("//")[:-1]
    for hmm in hmms:
        namepart = hmm.split("NAME  ")[1]
        name = namepart.split("\n")[0]
        lengthpart = hmm.split("LENG  ")[1]
        length = lengthpart.split("\n")[0]
        lengths[name] = int(length)
    return lengths

# DEAD FUNCTIONS
# these only exist so that the mapping to new functions is easier to do
def get_all_features_of_type(_seq_record, _types):
    raise RuntimeError("get_all_features_of_type(record, types) called, did you mean record.get_*()")

def get_cds_features_within_clusters(_seq_record):
    raise RuntimeError("get_withincluster_cds_features(record) called, use record.get_cds_features_within_clusters()")

def get_withincluster_cds_features(_seq_record):
    raise RuntimeError("get_withincluster_cds_features(record) called, use record.get_cds_features_within_clusters()")

def get_gene_id(_feature):
    raise RuntimeError("using get_gene_id(feature), did you mean feature.get_name() or feature.unique_id")

def get_aa_translation(_seq_record, _feature):
    raise RuntimeError("get_aa_translation(record, feature) called, use record.get_aa_translation(feature)")

def get_cluster_type(_cluster):
    raise RuntimeError("get_cluster_type(cluster) called, did you mean cluster.get_product_string() or cluster.products?")

def get_cluster_by_nr(_seq_record, _queryclusternr):
    raise RuntimeError("get_cluster_type(seq_record, cluster_num) called, did you mean seq_record.get_cluster(cluster_num)")

def sort_features(_seq_record):
    raise RuntimeError("utils.sort_features(seq_record) called, did you mean sorted(seq_record.get_all_features())?")

def get_cluster_cds_features(_cluster, _seq_record):
    raise RuntimeError("utils.get_cluster_cds_features(cluster) called, did you mean cluster.cds_children?")

def get_aa_sequence(feature, to_stop=False):
    raise RuntimeError("get_aa_sequence(cds) called, did you mean cds.get_aa_sequence()?")

