# License: GNU Affero General Public License v3 or later
# A copy of GNU AGPL v3 should have been included in this software package in LICENSE.txt.

""" Thiopeptides detection module

"""

import logging
from typing import Optional

from antismash.common import path
from antismash.config.args import ModuleArgs

from .specific_analysis import specific_analysis, ThioResults
from .html_output import generate_details_div, generate_sidepanel, will_handle

NAME = "thiopeptides"
SHORT_DESCRIPTION = NAME.capitalize()


def check_options(_options):
    """ No options to check for conflicts """
    return []


def is_enabled(options):
    """ Only disabled if minimal and not re-enabled """
    return options.thiopeptides_enabled or not options.minimal


def get_arguments():
    """ Build arguments, no explicit options as enabled by default """
    args = ModuleArgs('Advanced options', 'thio', enabled_by_default=True)
    return args


def check_prereqs():
    """ Check prereqs
            hmmpfam2: used to find extra HMM hits not in hmm_detection
    """
    failure_messages = []
    for binary_name in ['hmmpfam2']:
        if path.locate_executable(binary_name) is None:
            failure_messages.append("Failed to locate executable for %r" %
                                    binary_name)
    return failure_messages


def regenerate_previous_results(results, record, _options) -> Optional[ThioResults]:
    """ Regenerate previous results from JSON format """
    if not results:
        return None
    results = ThioResults.from_json(results, record)
    logging.debug("Reusing Thiopeptide results: %d clusters contained %d total motifs",
                  len(results.clusters_with_motifs), len(results.motifs))
    return results


def run_on_record(record, previous_results, _options) -> ThioResults:
    """ Run the thiopeptide analysis on the given record"""
    if previous_results:
        return previous_results
    return specific_analysis(record)
