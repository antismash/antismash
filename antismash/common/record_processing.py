# License: GNU Affero General Public License v3 or later
# A copy of GNU AGPL v3 should have been included in this software package in LICENSE.txt.

import functools
import logging
import re
import os
from typing import List

import Bio
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
from helperlibs.bio import seqio

from antismash.common import gff_parser
from antismash.common.secmet import Record
from antismash.config import get_config, update_config

from .subprocessing import parallel_function
from .utils import generate_unique_id


def parse_input_sequence(filename, minimum_length=-1, start=-1, end=-1) -> List[Record]:
    """ Parse input records contained in a file

        Arguments:
            filename: the path of the file to read
            minimum_length: records with length less than this will be ignored
                            if not positive, all records are included

        Returns:
            A list of secmet.Record instances, one for each record in the file
    """
    logging.info('Parsing input sequence %r', filename)
    if not isinstance(minimum_length, int):
        raise TypeError("minimum_length must be an int")

    records = []
    if not os.path.exists(filename):
        msg = "Sequence file not found: %r" % filename
        logging.error(msg)
        raise ValueError(msg)

    try:
        record_list = list(seqio.parse(filename))
        if not record_list:
            raise RuntimeError('No records could be read from file %r' % filename)
        for record in record_list:
            if minimum_length < 1 \
                    or len(record.seq) >= minimum_length \
                    or 'contig' in record.annotations \
                    or 'wgs_scafld' in record.annotations \
                    or 'wgs' in record.annotations:
                records.append(record)
    except (ValueError, AssertionError) as err:
        logging.error('Parsing %r failed: %s', filename, err)
        raise
    except Exception as err:
        logging.error('Parsing %r failed with unhandled exception: %s',
                      filename, err)
        raise

    # before conversion to secmet records, trim if required
    if start > -1 or end > -1:
        if len(records) > 1:
            raise ValueError("--start and --end options cannot be used with multiple records")
        trim_sequence(records[0], start, end)
    return [Record.from_biopython(record) for record in records]


def check_content(sequence: Record) -> Record:
    """ Checks if the sequence of a record is correct for the input type. If not
        the record's skip flag will be marked.

        Arguments:
            record: the Record instance to check

        Returns:
            the Record instance provided
    """
    cdsfeatures = sequence.get_cds_features()
    cdsfeatures_with_translations = len([cds for cds in cdsfeatures if cds.translation])
    assert cdsfeatures_with_translations == len(cdsfeatures)
    if not isinstance(sequence.seq.alphabet, Bio.Alphabet.NucleotideAlphabet)\
            and not is_nucl_seq(sequence.seq):
        logging.error("Record %s is a protein record, skipping.", sequence.id)
        sequence.skip = "protein record"
    else:
        sequence.seq.alphabet = Bio.Alphabet.generic_dna
    return sequence


def ensure_gene_info(single_entry: bool, genefinding, sequence: Record) -> Record:
    """ Ensures the given record has CDS features with unique locus tags.
        CDS features are retrieved from GFF file or via genefinding, depending
        on antismash options.

        Records without CDS features will have their skip flag marked.

        Arguments:
            single_entry: whether gff_parser can ignore mismatching record ids
                          provided there's only one record provided here and in
                          the GFF file
            genefinding: the relevant run_on_record(record, options) function to
                         use for finding genes if no GFF file being used
            record: the Record instance to ensure CDS features for

        Returns:
            the Record instance provided
    """
    options = get_config()
    if sequence.skip:
        return sequence
    if not sequence.get_cds_features():
        if options.genefinding_gff3:
            logging.info("No CDS features found in record %r but GFF3 file provided, running GFF parser.", sequence.id)
            gff_parser.run(sequence, single_entry, options)
            if not sequence.get_cds_features():
                logging.error("Record %s has no genes even after running GFF parser, skipping.", sequence.id)
                sequence.skip = "No genes found"
                return sequence
        elif options.genefinding_tool != "none":
            logging.info("No CDS features found in record %r, running gene finding.", sequence.id)
            genefinding(sequence, options)
        if not sequence.get_cds_features():
            logging.info("No genes found, skipping record")
            sequence.skip = "No genes found"
            return sequence
    fix_locus_tags(sequence)
    return sequence


def pre_process_sequences(sequences, options, genefinding) -> List[Record]:
    """ hmm

        - gaps removed
        - record ids adjusted to be unique
        - record ids are valid

        Note: Record instances will be altered in-place.

        Arguments:
            sequences: the secmet.Record instances to process
            options: an antismash Config instance
            genefinding: the module to use for genefinding, must have
                         run_on_record() implemented

        Returns:
            A list of altered secmet.Record
    """
    logging.debug("Preprocessing %d sequences", len(sequences))
    # keep count of how many records matched filter
    matching_filter = 0

    for i, seq in enumerate(sequences):
        seq.record_index = i

    # keep sequences as clean as possible
    logging.debug("Sanitising record sequences")
    sequences = parallel_function(sanitise_sequence, ([record] for record in sequences))
    for record in sequences:
        if record.skip or not record.seq:
            logging.warning("Record %s has no sequence, skipping.", record.id)

    # run the filter
    for sequence in sequences:
        if options.limit_to_record and options.limit_to_record != sequence.id:
            sequence.skip = "did not match filter: %s" % options.limit_to_record
        else:
            matching_filter += 1

    # Check if seq_records have appropriate content
    sequences = parallel_function(check_content, ([sequence] for sequence in sequences))

    if options.limit_to_record:
        limit = options.limit_to_record
        if matching_filter == 0:
            logging.error("No sequences matched filter: %s", limit)
            raise ValueError("No sequences matched filter: %s" % limit)
        elif matching_filter != len(sequences):
            logging.info("Skipped %d sequences not matching filter: %s",
                         len(sequences) - matching_filter, limit)

    # catch WGS master or supercontig entries
    if records_contain_shotgun_scaffolds(sequences):
        raise RuntimeError("Incomplete whole genome shotgun records are not supported")

    # Now remove small contigs < minimum length again
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
                    warned = True
                sequence.skip = "skipping all but first {0} meaningful records (--limit {0}) ".format(options.limit)

    options = update_config({"triggered_limit": warned})  # TODO is there a better way

    # Check if no duplicate locus tags / gene IDs are found
    ensure_no_duplicate_gene_ids(sequences)

    # Check GFF suitability
    single_entry = False
    if options.genefinding_gff3:
        single_entry = gff_parser.check_gff_suitability(options, sequences)

    all_record_ids = {seq.id for seq in sequences}
    # Ensure all records have unique names
    if len(all_record_ids) < len(sequences):
        all_record_ids = set()
        for record in sequences:
            if record.id in all_record_ids:
                record.original_id = record.id
                record.id = generate_unique_id(record.id, all_record_ids)[0]
            all_record_ids.add(record.id)
        assert len(all_record_ids) == len(sequences), "%d != %d" % (len(all_record_ids), len(sequences))
    # Ensure all records have valid names
    for record in sequences:
        fix_record_name_id(record, all_record_ids)

    partial = functools.partial(ensure_gene_info, single_entry, genefinding.run_on_record)
    sequences = parallel_function(partial, ([sequence] for sequence in sequences))

    return sequences


def sanitise_sequence(record) -> Record:
    """ Ensures all sequences use N for gaps instead of -, and that all other
        characters are A, C, G, T, or N

        Arguments:
            records: the secmet.Records to alter

        Returns:
            None
    """
    has_real_content = False
    sanitised = []
    for char in record.seq.upper():
        if char == "-":
            continue
        elif char in "ACGT":
            sanitised.append(char)
            has_real_content = True
        else:
            sanitised.append("N")
    record.seq = Seq("".join(sanitised), alphabet=record.seq.alphabet)
    if not has_real_content:
        record.skip = "contains no sequence"
    return record


def trim_sequence(record, start, end) -> SeqRecord:
    """ Trims a record to the range given

        Arguments:
            record: the Bio.SeqRecord to trim
            start: the start position (inclusive)
            end: the end position (exclusive)

        Returns:
            A new, shortened Bio.SeqRecord instance
    """
    if start >= len(record):
        raise ValueError('Specified analysis start point of %r is outside record' % start)
    if end > len(record):
        raise ValueError('Specified analysis end point of %r is outside record' % end)
    if end > -1 and end <= start:
        raise ValueError("Trim region start cannot be greater than or equal to end")

    if start < 0:
        start = 0
    if end < 0:
        end = len(record)
    return record[start:end]


def is_nucl_seq(sequence) -> bool:
    """ Determines if a sequence is a nucleotide sequence based on content.

        Arguments:
            sequence: the sequence to check, either a string or Bio.Seq

        Returns:
            True if less than 20% of bases are not a,c,g,t or n
    """
    other = str(sequence).lower()
    for char in "acgtn":
        other = other.replace(char, "")
    return len(other) < 0.2 * len(sequence)


def records_contain_shotgun_scaffolds(records) -> bool:
    """ Check if given records contain a WGS master record or supercontig record

        Arguments:
            records: an iterable of secmet.Record

        Returns:
            True if one of the given records is a WGS or supercontig record
    """
    for record in records:
        if 'wgs_scafld' in record.annotations \
                or 'wgs' in record.annotations \
                or 'contig' in record.annotations:
            return True
    return False


def ensure_no_duplicate_gene_ids(sequences) -> None:
    """ Ensures that every CDS across all sequences has a unique id

        Arguments:
            sequences: the secmet.Record instances to process

        Returns:
            None
    """
    all_ids = set()
    for sequence in sequences:
        for cdsfeature in sequence.get_cds_features():
            name = cdsfeature.get_name()
            if name in all_ids:
                name, _ = generate_unique_id(name[:8], all_ids, start=1)
            if cdsfeature.product is None:
                cdsfeature.product = name
            cdsfeature.locus_tag = name
            all_ids.add(name)


def fix_record_name_id(record, all_record_ids) -> None:
    """ Changes a record's name and id to be no more than 16 characters long,
        so it can be used in GenBank files.

        If record name is too long, the prefix c000X is used

        Arguments:
            record: the record to alter
            all_record_ids: a set of all known record ids

        Returns:
            None
    """

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
            contig_no = record.record_index

        return "c{ctg:05d}_{origid}..".format(ctg=contig_no, origid=idstring[:7])

    old_id = record.id

    if len(record.id) > 16:
        # Check if it is a RefSeq accession number like NZ_AMZN01000079.1 that
        # is too long just because of the version number behind the dot
        if (record.id[-2] == "." and
                record.id.count(".") == 1 and
                len(record.id.partition(".")[0]) <= 16 and
                record.id.partition(".")[0] not in all_record_ids):
            record.id = record.id.partition(".")[0]
            all_record_ids.add(record.id)
        else:  # Check if the ID suggested by _shorten_ids is unique
            if _shorten_ids(old_id) not in all_record_ids:
                name = _shorten_ids(old_id)
            else:
                name, _ = generate_unique_id(record.id[:12], all_record_ids, max_length=16)
            record.id = name
            all_record_ids.add(name)

        logging.warning('Fasta header too long: renamed "%s" to "%s"', old_id, record.id)

    if len(record.name) > 16:
        record.name = _shorten_ids(record.name)

    if 'accession' in record.annotations and \
       len(record.annotations['accession']) > 16:
        acc = record.annotations['accession']

        record.annotations['accession'] = _shorten_ids(acc)

    # Remove illegal characters from name: otherwise, file cannot be written
    illegal_chars = set('''!"#$%&()*+,:;=>?@[]^`'{|}/ ''')
    for char in record.id:
        if char in illegal_chars:
            record.id = record.id.replace(char, "")
    for char in record.name:
        if char in illegal_chars:
            record.name = record.name.replace(char, "")

    if not record.original_id and old_id != record.id:
        record.original_id = old_id


def fix_locus_tags(seq_record) -> None:  # TODO should be part of secmet
    "Fix CDS feature that don't have a locus_tag, gene name or protein id"
    for feature in seq_record.get_cds_features():
        # Fix locus tags, gene names or protein IDs if they contain illegal chars
        illegal_chars = set('''!"#$%&()*+,:; \r\n\t=>?@[]^`'{|}/ ''')
        for attr in ["locus_tag", "gene", "protein_id"]:
            val = getattr(feature, attr)
            if not val or not set(val).intersection(illegal_chars):
                continue
            for char in val:
                if char in illegal_chars:
                    val = val.replace(char, "_")
            logging.critical("%s altered in fix_locus_tags: %s->%s", attr, getattr(feature, attr), val)
            setattr(feature, attr, val)
