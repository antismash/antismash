# License: GNU Affero General Public License v3 or later
# A copy of GNU AGPL v3 should have been included in this software package in LICENSE.txt.

"""NRPS/PKS analysis module

"""
import datetime
from glob import glob
import os
import logging

from antismash.common import path, subprocessing
from antismash.config.args import ModuleArgs

from .html_output import will_handle, generate_sidepanel, generate_details_div, generate_js_domains
from .results import NRPS_PKS_Results
from .specific_analysis import specific_analysis

NAME = "nrps_pks"

SHORT_DESCRIPTION = "NRPS/PKS analysis"

# The tuple is the name of the binary and whether it is an optional requirement
_required_binaries = [
    ('hmmsearch', False),
    ('muscle', False),
    ('java', False),
]

def check_prereqs():
    failure_messages = []
    for binary_name, optional in _required_binaries:
        if path.locate_executable(binary_name) is None and not optional:
            failure_messages.append("Failed to locate executable for %r" %
                                    binary_name)

    return failure_messages


def get_arguments():
    args = ModuleArgs('NRPS/PKS options', 'np', enabled_by_default=True)
    args.add_analysis_toggle('--np-dummy',
                             dest='dummy',
                             action='store_true',
                             default=False,
                             help="placeholder")
    return args


def check_options(options):
    return []


def regenerate_previous_results(json, record, options):
    return NRPS_PKS_Results.from_json(json)


def is_enabled(options):
    return not options.minimal or options.nrps_pks_enabled


def run_on_record(record, results, options):
    if not results:
        results = NRPS_PKS_Results(record.id)
        specific_analysis(record, results, options)
    return results
