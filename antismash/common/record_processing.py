# License: GNU Affero General Public License v3 or later
# A copy of GNU AGPL v3 should have been included in this software package in LICENSE.txt.

""" A collection of functions for reading, processing, and sanitising records.
"""

import functools
import logging
import os
import re
from typing import Any, Dict, Callable, Iterable, List, Set, Tuple
import warnings

from Bio.Seq import Seq, UndefinedSequenceError
from Bio.SeqRecord import SeqRecord
from helperlibs.bio import seqio

from antismash.common import gff_parser
from antismash.common.errors import AntismashInputError
from antismash.common.secmet import Record
from antismash.common.secmet.errors import SecmetInvalidInputError
from antismash.config import get_config, update_config, Config, ConfigType
from antismash.custom_typing import AntismashModule

from .subprocessing import parallel_function


_INVALID_ID_CHARS = """!"#$%&()*+,:;=>?@[]^`'{|}/ """


def _strict_parse(filename: str) -> List[SeqRecord]:
    """ Parses the input record with extra wrappers to catch biopython warnings
        as errors.

        Arguments:
            filename: the name of the file to parse

        Returns:
            a list of SeqRecords parsed
    """
    filter_messages = [
        r".*invalid location.*",
        r".*Expected sequence length.*",
        r".*Couldn't parse feature location.*",
        r".*double-quote characters.*should be escaped.*",
    ]
    try:
        # prepend warning filters to raise exceptions on certain messages
        for message in filter_messages:
            warnings.filterwarnings("error", message=message)
        records = list(seqio.parse(filename))
    except Exception as err:
        message = str(err)
        # strip the "Ignoring" part, since it's not being ignored
        if message.startswith("Ignoring invalid location"):
            message = message[9:]
        if isinstance(err, AttributeError):
            message = f"error within parsing library: {message}"
        logging.error('Parsing %r failed: %s', filename, message)
        raise AntismashInputError(message) from err
    finally:
        # remove the new warning filters (functions in at least 3.5 and 3.6)
        warnings.filters = warnings.filters[len(filter_messages):]

    if not records:
        raise AntismashInputError(f"no valid records found in file {filename}")
    return records


def _strip_invalid_chars_from_ids(records: Iterable[Record]) -> None:
    def strip(value: str) -> str:
        return "".join(v for v in value if v not in _INVALID_ID_CHARS)

    for record in records:
        stripped_id = strip(record.id)
        # store the original name, if it hasn't already been modified
        if stripped_id != record.id and not record.original_id:
            record.original_id = record.id
        record.id = stripped_id
        record.name = strip(record.name)


def parse_input_sequence(filename: str, taxon: str = "bacteria", minimum_length: int = -1,
                         start: int = -1, end: int = -1, gff_file: str = "",
                         ignore_invalid_records: bool = False) -> List[Record]:
    """ Parse input records contained in a file

        Arguments:
            filename: the path of the file to read
            taxon: the taxon of the input, e.g. 'bacteria', 'fungi'
            minimum_length: records with length less than this will be ignored
                            if not positive, all records are included
            start: a start location for trimming the sequence, or -1 to use all
            end: an end location for trimming the sequence, or -1 to use all
            gff_file: a GFF file to use for gene/CDS annotations
            ignore_invalid_records: whether to continue past invalid records or not

        Returns:
            A list of secmet.Record instances, one for each record in the file
    """
    logging.info('Parsing input sequence %r', filename)
    if not isinstance(minimum_length, int):
        raise TypeError("minimum_length must be an int")

    records: List[SeqRecord] = []

    for record in _strict_parse(filename):
        # in no case allow record identifiers to be longer than the filesystem's
        # maximum file length, allowing for region numbers and extension
        if len(record.id) > os.pathconf("/", "PC_NAME_MAX") - len(".region000.gbk"):
            raise AntismashInputError(f"record identifier too long for file system: {record.id}")

        if records_contain_shotgun_scaffolds([record]):
            raise AntismashInputError("incomplete whole genome shotgun records are not supported")

        try:
            record.seq[0]
        except (IndexError, UndefinedSequenceError):
            raise AntismashInputError(f"record contains no sequence information: {record.id}")

        if minimum_length < 1 or len(record.seq) >= minimum_length:
            records.append(record)

    # if no records are left, that's a problem
    if not records:
        raise AntismashInputError(f"all input records smaller than minimum length ({minimum_length})")

    for record in records:
        if not Record.is_nucleotide_sequence(record.seq):
            raise AntismashInputError(f"protein records are not supported: {record.id}")

    # before conversion to secmet records, trim if required
    if start > -1 or end > -1:
        if len(records) > 1:
            raise AntismashInputError("--start and --end options cannot be used with multiple records")
        if start != -1 and end != -1 and start > end and records[0].annotations.get("topology") == "circular":
            raise AntismashInputError(
                "--start and --end cannot be used for "
                "a cross-origin section of a circular record"
            )
        records[0] = trim_sequence(records[0], max(start, 0), min(len(records[0]), end))

    # add GFF features before conversion, if relevant
    if gff_file:
        logging.debug("Loading annotations from GFF file")
        # check GFF suitability first
        try:
            gff_parser.check_gff_suitability(gff_file, records)
        except AntismashInputError:
            raise
        except Exception as err:
            # avoid swallowing details if possible
            if str(err):
                logging.error(err)
            raise AntismashInputError("could not parse records from GFF3 file") from err
        gff_parser.update_records(gff_file, records)

    # remove any previous or obselete antiSMASH annotations to minimise incompatabilities
    for record in records:
        strip_record(record)

    logging.debug("Converting records from biopython to secmet")
    converted_records = []
    for record in records:
        try:
            converted_records.append(Record.from_biopython(record, taxon, discard_antismash_features=True))
        except SecmetInvalidInputError as err:
            if ignore_invalid_records:
                logging.warning("Ignoring invalid record %r %s", record.id, err)
            else:
                raise AntismashInputError(str(err)) from err
    if not converted_records:
        raise AntismashInputError(f"no valid records found in file {filename}")

    # if parsable by secmet, it has a better context on what to strip, so run
    # the secmet stripping to ensure there's no surprises
    for record in converted_records:
        record.strip_antismash_annotations()

    return converted_records


def strip_record(record: SeqRecord) -> SeqRecord:
    """ Discard antismash specific features and feature qualifiers

        Arguments:
            record: the SeqRecord object to clean

        Returns
            the same SeqRecord instance that was given
    """
    logging.debug("Stripping antiSMASH features and annotations from record: %s", record.id)
    old_as_types = {"cluster", "cluster_border", "cluster_core"}
    current_as_types = {"protocluster", "proto_core", "cand_cluster", "subregion",
                        "region", "PFAM_domain", "aSDomain", "promoter"}

    as_cds_qualifiers = {
        "aSProdPred",  # <= aS4
        "sec_met",  # <= aS4, aS5-alpha
        "gene_functions",  # aS5-beta
        "sec_met_domains",  # aS5-beta
        "gene_kind",  # aS5-beta
        "NRPS_PKS",  # as5.0.0
    }

    # remove any record-level comment
    if "structured_comment" in record.annotations:
        record.annotations["structured_comment"].pop("antiSMASH-Data", None)
        if not record.annotations["structured_comment"]:
            record.annotations.pop("structured_comment")

    # remove relevant features and qualifiers
    kept_features = []
    for feature in record.features:
        # remove anything explicitly added by antismash, this will catch shared types (e.g. CDS_motif)
        if "antismash" in [tool.lower() for tool in feature.qualifiers.get("tool", [])]:
            continue
        if feature.qualifiers.get("aSTool"):
            continue
        # then check types explicitly as they can be too old to contain the tool
        if feature.type in old_as_types:
            continue  # remove to avoid confusion
        if feature.type in current_as_types:
            continue

        # old CDS_motif annotations didn't include a "tool" qualifier, so clean them up
        if feature.type == "CDS_motif":
            # remove NRPS/PKS domain annotations
            notes = [note for note in feature.qualifiers.get("note", []) if not note.startswith("NRPS/PKS Motif")]
            # ensure notes isn't left in but empty
            if not notes:
                feature.qualifiers.pop("note", [])
            else:
                feature.qualifiers["note"] = notes
            # remove any feature that now has no qualifiers
            if not feature.qualifiers:
                continue

        # since CDS features can have a lot of old annotations, clean them up too
        if feature.type == "CDS":
            for qual in as_cds_qualifiers:
                feature.qualifiers.pop(qual, None)
            notes = feature.qualifiers.get("note")
            if notes:  # smcogs
                notes = [note for note in notes if not note.startswith("smCOG")]
                feature.qualifiers["note"] = notes

        kept_features.append(feature)

    record.features = kept_features

    return record


def ensure_cds_info(genefinding: Callable[[Record, Any], None], sequence: Record,
                    **kwargs: Dict[str, Any]) -> Record:
    """ Ensures the given record has CDS features with unique locus tags.
        CDS features are retrieved from GFF file or via genefinding, depending
        on antismash options.

        Records without CDS features will have their skip flag marked.

        Arguments:
            genefinding: the relevant run_on_record(record, options) function to
                         use for finding genes if no GFF file being used
            record: the Record instance to ensure CDS features for

        Returns:
            the Record instance provided
    """
    if sequence.skip:
        return sequence
    options = get_config(no_defaults=True)
    if len(options) == 0:  # inside a parallel function where config doesn't pickle correctly
        new = Config(kwargs)
        assert isinstance(new, ConfigType)
        options = new
    if not sequence.get_cds_features():
        if not options.genefinding_gff3 and options.genefinding_tool != "none":
            logging.info("No CDS features found in record %r, running gene finding.", sequence.id)
            try:
                genefinding(sequence, options)
            except ValueError as err:
                raise AntismashInputError(str(err)) from err
        if not sequence.get_cds_features():
            logging.info("No genes found, skipping record %r", sequence.id)
            sequence.skip = "No genes found"
            return sequence
    return sequence


def filter_records_by_name(sequences: List[Record], target: str) -> None:
    """ Mark records as skipped if their id does not match the given target or
        they are above .

        If the target is an empty string, all records will match.
        If not records match, an error will be raised.

        Arguments:
            sequences: the Records to filter
            target: the name to match, must be exact

        Returns:
            None
    """
    if not target:
        return

    logging.debug("Limiting to record id: %s", target)

    # run the filter
    matching_filter = 0
    for sequence in sequences:
        if sequence.id != target:
            sequence.skip = f"did not match filter: {target}"
        else:
            matching_filter += 1

    if matching_filter == 0:
        logging.error("No sequences matched filter: %s", target)
        raise AntismashInputError(f"no sequences matched filter: {target}")

    logging.info("Skipped %d sequences not matching filter: %s",
                 len(sequences) - matching_filter, target)


def filter_records_by_count(records: List[Record], maximum: int) -> bool:
    """ Mark all records after the largest 'maximum' non-skipped records as skipped.

        If maximum is -1, no records will be skipped due to count.

    Arguments:
        records: the Records to filter
        maximum: the maximum number of records to run

    Returns:
        True if any records were marked as skipped due to hitting the limit
    """
    if maximum == -1 or maximum > len(records):
        return False

    # sort by decreasing length, breaking ties by keeping existing ordering
    records = [pair[1] for pair in sorted(enumerate(records), key=lambda x: (-len(x[1]), x[0]))]

    limit_hit = False
    meaningful = 0
    for record in records:
        if record.skip:  # don't count any records that are already skipped
            continue
        meaningful += 1
        if meaningful > maximum:
            limit_hit = True
            record.skip = f"skipping all but largest {maximum} meaningful records (--limit) "

    return limit_hit


def pre_process_sequences(sequences: List[Record], options: ConfigType, genefinding: AntismashModule) -> List[Record]:
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

    # catch WGS master or supercontig entries
    if records_contain_shotgun_scaffolds(sequences):
        raise AntismashInputError("incomplete whole genome shotgun records are not supported")

    for i, seq in enumerate(sequences):
        seq.record_index = i + 1  # 1-indexed

    checking_required = not (options.reuse_results or options.skip_sanitisation)

    # keep sequences as clean as possible and make sure they're valid
    if checking_required:
        logging.debug("Sanitising record ids and sequences")
        _strip_invalid_chars_from_ids(sequences)
        # Ensure all records have unique names
        all_record_ids = {seq.id for seq in sequences}
        if len(all_record_ids) < len(sequences):
            all_record_ids = set()
            for record in sequences:
                if record.id in all_record_ids:
                    if not record.original_id:
                        record.original_id = record.id
                    record.id = generate_unique_id(record.id, all_record_ids)[0]
                all_record_ids.add(record.id)
            assert len(all_record_ids) == len(sequences), f"{len(all_record_ids)} != {len(sequences)}"
        # Ensure all records have valid names
        for record in sequences:
            fix_record_name_id(record, all_record_ids, options.allow_long_headers)
        if len(sequences) == 1:
            sequences = [sanitise_sequence(sequences[0])]
        else:
            sequences = parallel_function(sanitise_sequence, ([record] for record in sequences))

    for record in sequences:
        if record.skip or not record.seq:
            logging.warning("Record %s has no sequence, skipping.", record.id)
        if not record.id:
            raise AntismashInputError("record has no name")

    # skip anything not matching the filter
    filter_records_by_name(sequences, options.limit_to_record)

    # Now remove small contigs < minimum length again
    logging.debug("Removing sequences smaller than %d bases", options.minlength)
    for sequence in sequences:
        if len(sequence.seq) < options.minlength:
            sequence.skip = f"smaller than minimum length ({options.minlength})"

    # Make sure we don't waste weeks of runtime on huge records, unless requested by the user
    limit_hit = filter_records_by_count(sequences, options.limit)
    if limit_hit:
        logging.warning("Only analysing the first %d records (increase via --limit)", options.limit)
    update_config({"triggered_limit": limit_hit})

    if checking_required:
        # ensure CDS features have all relevant information
        logging.debug("Ensuring records have CDS features with all required information")
        assert hasattr(genefinding, "run_on_record")
        genefinding_opts = {key: val for key, val in options if key.startswith("genefinding")}
        genefinding_opts["taxon"] = options.taxon
        partial = functools.partial(ensure_cds_info, genefinding.run_on_record, **genefinding_opts)
        sequences = parallel_function(partial, ([sequence] for sequence in sequences))

    if all(sequence.skip for sequence in sequences):
        raise AntismashInputError("all records skipped")

    return sequences


def sanitise_sequence(record: Record) -> Record:
    """ Ensures all sequences use N for gaps instead of -, and that all other
        characters are A, C, G, T, or N

        Arguments:
            records: the secmet.Records to alter

        Returns:
            the same Record instance as given
    """
    has_real_content = False
    sanitised = []
    for char in record.seq.upper():
        if char == "-":
            continue
        if char in "ACGT":
            sanitised.append(char)
            has_real_content = True
        else:
            sanitised.append("N")
    record.seq = Seq("".join(sanitised))
    if not has_real_content:
        record.skip = "contains no sequence"
    return record


def trim_sequence(record: SeqRecord, start: int, end: int) -> SeqRecord:
    """ Trims a record to the range given

        Arguments:
            record: the Bio.SeqRecord to trim
            start: the start position (inclusive)
            end: the end position (exclusive)

        Returns:
            A new, shortened Bio.SeqRecord instance
    """
    if start >= len(record):
        raise ValueError(f"Specified analysis start point of {start!r} is outside record")
    if end > len(record):
        raise ValueError(f"Specified analysis end point of {end!r} is outside record")
    if -1 < end <= start:
        raise ValueError("Trim region start cannot be greater than or equal to end")

    start = max(start, 0)
    if end < 0:
        end = len(record)
    new = record[start:end]
    new.annotations = record.annotations  # preserve old annotations that biopython ignores
    return new


def records_contain_shotgun_scaffolds(records: List[Record]) -> bool:
    """ Check if given records contain a WGS master record or supercontig record

        Arguments:
            records: an iterable of secmet.Record

        Returns:
            True if one of the given records is a WGS or supercontig record
    """
    for record in records:
        defined = True
        try:
            record.seq[0]
        except UndefinedSequenceError:
            defined = False
        except IndexError:
            raise AntismashInputError(f"record contains no sequence information: {record.id}")
        if not defined and any(key in record.annotations for key in [
            "wgs_scafld",
            "wgs",
            "contig",
        ]):
            return True
    return False


def fix_record_name_id(record: Record, all_record_ids: Set[str],
                       allow_long_names: bool = False) -> None:
    """ Changes a record's name and id to be no more than 16 characters long,
        so it can be used in GenBank files.

        If record name is too long, the prefix c000X is used

        Arguments:
            record: the record to alter
            all_record_ids: a set of all known record ids
            allow_long_names: whether names longer than 16 characters are allowed

        Returns:
            None
    """

    def _shorten_ids(idstring: str) -> str:
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
            assert isinstance(record.record_index, int)
            contig_no = record.record_index

        return f"c{contig_no:05d}_{idstring[:7]}.."

    old_id = record.id

    if len(record.id) > 16 and not allow_long_names:
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

    if len(record.name) > 16 and not allow_long_names:
        record.name = _shorten_ids(record.name)

    if 'accession' in record.annotations and \
       len(record.annotations['accession']) > 16:
        acc = record.annotations['accession']

        record.annotations['accession'] = _shorten_ids(acc)

    if not record.original_id and old_id != record.id:
        record.original_id = old_id


def generate_unique_id(prefix: str, existing_ids: Set[str], start: int = 0,
                       max_length: int = -1) -> Tuple[str, int]:
    """ Generate a identifier of the form prefix_num, e.g. seq_15.

        Does *not* add the generated prefix to current identifiers

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
    max_length = int(max_length)

    name = f"{prefix}_{counter}"
    while name in existing_ids:
        counter += 1
        name = f"{prefix}_{counter}"
    if 0 < max_length < len(name):
        raise RuntimeError(f"Could not generate unique id for {prefix} after {counter - start} iterations")
    return name, counter
