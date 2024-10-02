# License: GNU Affero General Public License v3 or later
# A copy of GNU AGPL v3 should have been included in this software package in LICENSE.txt.

""" A module to find genes in a record using external tools """

import logging
from typing import List

from antismash.common.errors import AntismashInputError
from antismash.common.secmet import Record
from antismash.config import ConfigType
from antismash.config.args import ModuleArgs, ReadableFullPathAction

from .run_prodigal import run_prodigal


NAME = "genefinding"
SHORT_DESCRIPTION = "Genefinding with Prodigal for prokaryotes"


def get_arguments() -> ModuleArgs:
    """ Construct args with options for which genefinding method to use (if any)
    """
    args = ModuleArgs('Gene finding options (ignored when ORFs are annotated)',
                      'genefinding', basic_help=True)
    args.add_option('tool',
                    dest='tool',
                    default='error',
                    choices=['prodigal', 'prodigal-m', 'none', 'error'],
                    type=str,
                    help="Specify algorithm used for gene finding: "
                         "Prodigal, Prodigal Metagenomic/Anonymous mode, or none."
                         " The 'error' option will raise an error if genefinding is attempted."
                         " The 'none' option will not run genefinding."
                         " (default: %(default)s).")
    args.add_option('gff3',
                    dest='gff3',
                    default="",
                    action=ReadableFullPathAction,
                    type=str,
                    metavar="GFF3_FILE",
                    help="Specify GFF3 file to extract features from.")
    return args


def check_prereqs(options: ConfigType) -> List[str]:
    """ Make sure the external tools to use are available """
    failure_messages: List[str] = []
    if options.genefinding_tool in ['none']:
        return failure_messages
    binaries: List[str] = []
    if options.check_prereqs_only:
        binaries = ["prodigal"]
    elif options.genefinding_tool in ['prodigal', 'prodigal-m']:
        binaries = ['prodigal']
    for binary_name in binaries:
        if binary_name not in options.executables:
            failure_messages.append(f"Failed to locate executable for {binary_name}")

    return failure_messages


def check_options(options: ConfigType) -> List[str]:
    """ Check that fungal sequences aren't using bacterial genefinding
        and vice versa.
    """
    errors = []
    if options.taxon == "fungi" and options.genefinding_tool not in ["none", "error"]:
        errors.append("Fungi taxon must provide gene annotations in GenBank or GFF3 format.")
    return errors


def is_enabled(options: ConfigType) -> bool:
    """ Considered enabled in the case of the tool being 'error', because it
        should throw an error in that case.
    """
    return options.genefinding_tool != "none"


def run_on_record(record: Record, options: ConfigType) -> None:
    """ Find genes in a Record using prodigal.
        Genes will be added to the record as they are found.
    """
    if options.genefinding_tool in ["none"]:
        return None

    if options.taxon == "fungi":
        raise AntismashInputError(
            f"Fungal record {record.id} contains no genes and "
            "antiSMASH does not provide a fungal gene finder."
        )

    if options.genefinding_tool == "error":
        raise AntismashInputError(f"Record {record.id} contains no genes and no genefinding tool specified")

    if options.genefinding_tool in ["prodigal", "prodigal-m"]:
        logging.debug("Running prodigal based genefinding")
        return run_prodigal(record)

    raise ValueError(f"Unknown genefinding tool: {options.genefinding_tool}")
