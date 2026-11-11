# License: GNU Affero General Public License v3 or later
# A copy of GNU AGPL v3 should have been included in this software package in LICENSE.txt.

"""
Genbank output module
"""

from argparse import BooleanOptionalAction
from datetime import datetime
import logging
import os

from Bio.SeqRecord import SeqRecord
from helperlibs.bio import seqio as hl_seqio

from antismash.common.module_results import ModuleResults
from antismash.common.secmet import Record, Region
from antismash.custom_typing import AntismashModule
from antismash.config import ConfigType
from antismash.config.args import ModuleArgs


NAME = "gbk"
SHORT_DESCRIPTION = "Genbank output"


def add_antismash_comments(records: list[tuple[Record, SeqRecord]], options: ConfigType) -> None:
    """ Add antismash meta-annotation to records for genbank output

        Arguments:
            records: a list of Record, SeqRecord pairs
            options: antismash options

        Returns:
            None

    """
    if not records:
        return
    base_comment = {
        "Version": str(options.version),
        "Run date": str(datetime.now().strftime('%Y-%m-%d %H:%M:%S'))
    }
    # include start/end details if relevant
    if options.start != -1 or options.end != -1:
        start = 1 if options.start == -1 else options.start
        # start/end is only valid for single records, as per record_processing
        assert len(records) == 1
        end = len(records[0][0].seq) if options.end == -1 else options.end
        base_comment.update({
            "NOTE": "This is an extract from the original record!",
            "Starting at": str(start),
            "Ending at": str(end),
        })

    for record, bio_record in records:
        comment = dict(base_comment)
        if record.original_id:
            comment["Original ID"] = record.original_id

        if "structured_comment" not in bio_record.annotations:
            bio_record.annotations["structured_comment"] = {}

        bio_record.annotations["structured_comment"][f"{options.branding}-Data"] = comment


def check_options(_options: ConfigType) -> list[str]:
    """ Check options are valid """
    return []


def check_prereqs(_options: ConfigType) -> list[str]:
    """ Check prerequisites """
    return []


def get_arguments() -> ModuleArgs:
    """ Builds the arguments for the HMTL output module """
    args = ModuleArgs("Output options", "gbk", enabled_by_default=True)
    args.add_option(
        "summary",
        dest="gbk_summary",
        action=BooleanOptionalAction,
        default=True,
        help=("Output an annotated genbank file containing all input records (default: %(default)s"),
    )
    args.add_option(
        "regions",
        dest="gbk_regions",
        action=BooleanOptionalAction,
        default=True,
        help=("Output an annotated genbank file per region (default: %(default)s"),
    )
    args.add_option(
        "compress-summary",
        dest="gbk_compress_summary",
        action=BooleanOptionalAction,
        default=False,
        help=("Output an annotated genbank file containing all input records (default: %(default)s"),
    )
    return args


def get_region_file_path(record: Record, region: Region, options: ConfigType) -> str | None:
    """ Builds the path for a region's genbank file, if it would be written.

        Arguments:
            record: the record the region belongs to
            region: the region for which to build a path
            options: the options to use when building the path

        Returns:
            the path of the region file, or None if outputs are disabled
    """
    if not (is_enabled(options) and options.gbk_regions):
        return None
    return f"{record.id}.region{region.get_region_number():03d}.gbk"


def get_summary_file_path(options: ConfigType) -> str | None:
    """ Builds the path for the summary file, if it would be written

        Arguments:
            options: the config to build from

        Returns:
            the path of the summary file, or None if outputs are disabled
    """
    if not (is_enabled(options) and options.gbk_summary):
        return None
    summary = os.path.join(options.output_dir, options.output_basename + ".gbk")
    if options.gbk_compress_summary:
        summary += ".gz"
    return summary


def is_enabled(options: ConfigType) -> bool:
    """ Returns true if the module is enabled """
    return all([
        not options.is_module_excluded(NAME),
        options.gbk_enabled,
        any([
            options.gbk_summary,
            options.gbk_regions,
        ]),
    ])


def prepare_data(_logging_only: bool = False) -> list[str]:
    """ Rebuild any dynamically buildable data """
    return []


def write(records: list[Record], _results: list[dict[str, ModuleResults]],
          options: ConfigType, _all_modules: list[AntismashModule]) -> None:
    """ Writes all results to file, where applicable.

        Arguments:
            records: the list of Records for which results exist
            results: a list of dictionaries containing all module results for records
            options: antismash config object
            all_modules: the modules used to generate the results

        Returns:
            None
    """
    bio_records = [record.to_biopython() for record in records]

    add_antismash_comments(list(zip(records, bio_records)), options)

    if options.gbk_regions:
        logging.debug("Writing genbank files for individual regions")
        for record, bio_record in zip(records, bio_records):
            for region in record.get_regions():
                filename = get_region_file_path(record, region, options)
                region.write_to_genbank(filename=filename, directory=options.output_dir,
                                        record=bio_record)

    # write records to an aggregate output
    if options.gbk_summary:
        summary_path = get_summary_file_path(options)
        logging.debug("Writing summary genbank file to '%s'", summary_path)
        hl_seqio.write(bio_records, summary_path)
