# License: GNU Affero General Public License v3 or later
# A copy of GNU AGPL v3 should have been included in this software package in LICENSE.txt.

""" A module to find genes in a record using external tools """

import logging
import os
from typing import List

from antismash.common.path import locate_executable
from antismash.common.secmet import Record
from antismash.config import get_config, ConfigType
from antismash.config.args import ModuleArgs

from .run_prodigal import run_prodigal
from .run_glimmerhmm import run_glimmerhmm


NAME = "genefinding"
SHORT_DESCRIPTION = "Genefinding with GlimmerHMM or Prodigal"


def get_arguments() -> ModuleArgs:
    """ Construct args with options for which genefinding method to use (if any)
    """
    args = ModuleArgs('Gene finding options (ignored when ORFs are annotated)', 'genefinding')
    args.add_option('tool',
                    dest='tool',
                    default='error',
                    choices=['glimmerhmm', 'prodigal', 'prodigal-m', 'none', 'error'],
                    type=str,
                    help="Specify algorithm used for gene finding: GlimmerHMM, "
                         "Prodigal, Prodigal Metagenomic/Anonymous mode, or none."
                         " The 'error' option will raise an error if genefinding is attempted."
                         " The 'none' option will not run genefinding."
                         " (default: %(default)s).")
    args.add_option('gff3',
                    dest='gff3',
                    default="",
                    type=str,
                    metavar="GFF3_FILE",
                    help="Specify GFF3 file to extract features from.")
    return args


def check_prereqs() -> List[str]:
    """ Make sure the external tools to use are available """
    failure_messages = []  # type: List[str]
    options = get_config()
    if options.genefinding_tool in ['none']:
        return failure_messages
    binaries = []  # type: List[str]
    if options.check_prereqs_only:
        binaries = ["prodigal", "glimmerhmm"]
    elif options.genefinding_tool in ['prodigal', 'prodigal-m']:
        binaries = ['prodigal']
    elif options.taxon == 'fungi':
        binaries = ['glimmerhmm']
    for binary_name in binaries:
        if not locate_executable(binary_name):
            failure_messages.append("Failed to locate executable for %r" % binary_name)

    return failure_messages


def check_options(options: ConfigType) -> List[str]:
    """ Check that fungal sequences aren't using bacterial genefinding
        and vice versa.
    """
    errors = []
    if options.genefinding_gff3:
        if not os.path.exists(options.genefinding_gff3):
            errors.append("Specified gff file does not exist: %s" % (
                    options.genefinding_gff3))
    if options.taxon == "fungi" and options.genefinding_tool not in ["glimmerhmm", "none", "error"]:
        errors.append("Fungi taxon must use glimmerhmm for genefinding if using genefinding")
    if options.taxon == "bacteria" and options.genefinding_tool == "glimmerhmm":
        errors.append("Bacteria taxon cannot use glimmerhmm for genefinding")
    return errors


def is_enabled(options: ConfigType) -> bool:
    """ Considered enabled in the case of the tool being 'error', because it
        should throw an error in that case.
    """
    return options.genefinding_tool != "none"


def run_on_record(record: Record, options: ConfigType) -> None:
    """ Find genes in a Record using glimmerhmm or prodigal.
        Genes will be added to the record as they are found.
    """
    if options.genefinding_tool == 'error':
        raise ValueError("Called find_genes, but genefinding disabled")

    if options.taxon == 'fungi':
        if options.genefinding_tool == ["none"]:
            return None
        assert options.genefinding_tool == "glimmerhmm"
        logging.debug("Running glimmerhmm genefinding")
        return run_glimmerhmm(record)
    elif options.genefinding_tool in ["prodigal", "prodigal-m"]:
        logging.debug("Running prodigal based genefinding")
        return run_prodigal(record, options)

    raise ValueError("Unknown genefinding tool: %s" % options.genefinding_tool)
