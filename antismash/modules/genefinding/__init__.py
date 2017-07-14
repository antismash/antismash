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
                      default='prodigal',
                      choices=['glimmer', 'prodigal', 'prodigal-m', 'none'],
                      type=str,
                      help="Specify algorithm used for gene finding: Glimmer, "
                           "Prodigal, Prodigal Metagenomic/Anonymous mode or none. (default: %(default)s).")
    args.add_argument('all-orfs',
                      dest='all_orfs',
                      action='store_true',
                      default=False,
                      help="Use all ORFs > 60 nucleotides instead of running genefinding or using a GFF.")
    args.add_argument('gff3',
                      dest='gff3',
                      default=False,
                      type=bool,
                      help="Specify GFF3 file to extract features from.")
    return args

def check_prereqs():
    failure_messages = []
    options = Config()
    if options.genefinding_tool == 'none':
        return failure_messages
    basedir = options.get('glimmer', {}).basedir
    for binary_name in ['long-orfs', 'extract', 'build-icm',
                                  'glimmer3', 'glimmerhmm', 'prodigal']:
        new_binary_name = os.path.join(basedir, binary_name)
        if not locate_executable(new_binary_name) and not locate_executable(binary_name):
            failure_messages.append("Failed to locate executable for %r" % binary_name)

    return failure_messages

def check_options(options):
    errors = []
    if options.genefinding_gff3:
        if not os.path.exists(options.genefinding_gff3):
            errors.append("Specified gff file does not exist: %s" % (
                    options.genefinding_gff3))
    return errors

def is_enabled(options):
    return options.genefinding_tool != "none" or options.genefinding_gff3
