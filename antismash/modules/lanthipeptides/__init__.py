# License: GNU Affero General Public License v3 or later
# A copy of GNU AGPL v3 should have been included in this software package in LICENSE.txt.

""" Lanthipeptides detection module

"""

import logging

import antismash.common.path as path
from antismash.config.args import ModuleArgs

from .config import get_config # TODO: what is this
from .specific_analysis import specific_analysis, LanthiResults
from .html_output import generate_details_div, generate_sidepanel, will_handle

NAME = "lanthipeptides"
SHORT_DESCRIPTION = NAME.capitalize()

def get_arguments():
    args = ModuleArgs('Advanced options', 'lanthi', enabled_by_default=True)
    return args

def check_options(options):
    return []

def is_enabled(options):
    return not (options.minimal and not options.lanthipeptides_enabled) # TODO: use --minimal and relevant --enable instead

def check_previous_results(results, record, options):
    if not results:
        return None
    results = LanthiResults.from_json(results, record)
    logging.debug("Reusing Lanthipeptide results: %d clusters contained %d total motifs",
                  len(results.clusters_with_motifs), len(results.motifs))
    return results

def check_prereqs():
    failure_messages = []
    for binary_name, optional in [('hmmpfam2', False), ('fimo', True)]:
        present = True
        if path.locate_executable(binary_name) is None:
            present = False
            if not optional:
                failure_messages.append("Failed to locate executable for %r" %
                                        binary_name)
        slot = '{}_present'.format(binary_name)
        conf = get_config()
        if hasattr(conf, slot):
            setattr(conf, slot, present)

    return failure_messages

def run_on_record(record, options):
    return specific_analysis(record)
