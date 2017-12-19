# License: GNU Affero General Public License v3 or later
# A copy of GNU AGPL v3 should have been included in this software package in LICENSE.txt.

"""Sactipeptides detection module

"""

from antismash.common import path
from antismash.config.args import ModuleArgs

from .specific_analysis import specific_analysis, SactiResults
from .html_output import generate_details_div, generate_sidepanel, will_handle

NAME = "sactipeptides"
SHORT_DESCRIPTION = "sactipeptide detection"


def check_prereqs():
    failure_messages = []
    for binary_name in ['hmmpfam2']:
        if not path.locate_executable(binary_name):
            failure_messages.append("Failed to locate executable for %r" %
                                    binary_name)
    return failure_messages


def get_arguments():
    """ Runs by default, but add minimal's --enable option """
    args = ModuleArgs('Advanced options', 'sacti', enabled_by_default=True)
    return args


def check_options(options):
    return []


def is_enabled(options):
    return not options.minimal or options.sactipeptides_enabled


def regenerate_previous_results(previous, record, options):
    return SactiResults.from_json(previous, record)


def run_on_record(record, results, options):
    return specific_analysis(record, options)
