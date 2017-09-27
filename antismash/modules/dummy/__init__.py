# License: GNU Affero General Public License v3 or later
# A copy of GNU AGPL v3 should have been included in this software package in LICENSE.txt.

"""reserves/creates some cmdline options that we need placeholders for"""

import argparse

from antismash.config.args import ModuleArgs

NAME = "dummy"
SHORT_DESCRIPTION = "options placeholders"

def get_arguments():
    args = ModuleArgs('Dummy options', '', override_safeties=True)
    args.add_option('--dummy-cassis',
                      dest='cassis',
                      action='store_true',
                      default=False,
                      help="Dummy only.")
    args.add_option('--dummy-borderpredict',
                      dest='borderpredict',
                      action='store_true',
                      default=False,
                      help="Dummy only.")
    args.add_option('--dummy-without-fimo',
                      dest='without_fimo',
                      action='store_true',
                      default=False,
                      help="Dummy only.")


    return args

def check_options(options):
    if options.without_fimo or options.borderpredict:
        raise ValueError("Dummy options can't be enabled")
    return []

def check_prereqs():
    """Check for prerequisites"""
    # No external dependencies
    return []

def is_enabled(options):
    return False

def regenerate_previous_results(previous, record, options):
    return None

def run_on_record(seq_record, results, options):
    raise NotImplementedError("Dummy module should never be run")
