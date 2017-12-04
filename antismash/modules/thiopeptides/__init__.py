# License: GNU Affero General Public License v3 or later
# A copy of GNU AGPL v3 should have been included in this software package in LICENSE.txt.

""" Thiopeptides detection module

"""

import logging

from antismash.common import path
from antismash.config.args import ModuleArgs

from .specific_analysis import specific_analysis, ThioResults
from .html_output import generate_details_div, generate_sidepanel, will_handle

NAME = "thiopeptides"
SHORT_DESCRIPTION = NAME.capitalize()


def check_options(options):
    return []


def is_enabled(options):
    return options.thiopeptides_enabled or not options.minimal


def get_arguments():
    args = ModuleArgs('Advanced options', 'thio', enabled_by_default=True)
    return args


def check_prereqs():
    failure_messages = []
    for binary_name, optional in [('hmmpfam2', False)]:
        if path.locate_executable(binary_name) is None and not optional:
            failure_messages.append("Failed to locate executable for %r" %
                                    binary_name)
    return failure_messages


def regenerate_previous_results(results, record, options):
    if not results:
        return None
    results = ThioResults.from_json(results, record)
    logging.debug("Reusing Thiopeptide results: %d clusters contained %d total motifs",
                  len(results.clusters_with_motifs), len(results.motifs))
    return results


def run_on_record(record, previous_results, options):
    if previous_results:
        return previous_results
    return specific_analysis(record, options)
