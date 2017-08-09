# License: GNU Affero General Public License v3 or later
# A copy of GNU AGPL v3 should have been included in this software package in LICENSE.txt.

import logging
import os

from antismash.config.args import Config
from antismash.common.path import locate_executable
from antismash.config.args import ModuleArgs
from antismash.modules.genefinding.genefinding import run_on_record

NAME = "genefinding"
SHORT_DESCRIPTION = NAME.capitalize()

def get_arguments():
    args = ModuleArgs('Gene finding options (ignored when ORFs are annotated)', 'genefinding')
    args.add_argument('tool',
                      dest='tool',
                      default='none',
                      choices=['glimmerhmm', 'prodigal', 'prodigal-m', 'all-orfs', 'none'],
                      type=str,
                      help="Specify algorithm used for gene finding: GlimmerHMM, "
                           "Prodigal, Prodigal Metagenomic/Anonymous mode, use"
                           " all ORFs > 60 nucleotides, or none."
                           " (default: %(default)s).")
    args.add_argument('gff3',
                      dest='gff3',
                      default="",
                      type=str,
                      help="Specify GFF3 file to extract features from.")
    return args

def check_prereqs():
    failure_messages = []
    options = Config()
    if options.genefinding_tool in ['none', 'all-orfs']:
        return failure_messages
    binaries = []
    if options.genefinding_tool in ['prodigal', 'prodigal-m']:
        binaries = ['prodigal']
    elif options.taxon == 'fungi':
        binaries = ['glimmerhmm']
    for binary_name in binaries:
        if not locate_executable(binary_name):
            failure_messages.append("Failed to locate executable for %r" % binary_name)

    return failure_messages

def check_options(options):
    errors = []
    if options.genefinding_gff3:
        if not os.path.exists(options.genefinding_gff3):
            errors.append("Specified gff file does not exist: %s" % (
                    options.genefinding_gff3))
    if options.taxon == "fungi" and options.genefinding_tool not in ["glimmerhmm", "none"]:
        errors.append("Fungi taxon must use glimmerhmm for genefinding if using genefinding")
    if options.taxon == "bacteria" and options.genefinding_tool == "glimmerhmm":
        errors.append("Bacteria taxon cannot use glimmerhmm for genefinding")
    return errors

def is_enabled(options):
    return options.genefinding_tool != "none" or options.genefinding_gff3
