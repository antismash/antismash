# License: GNU Affero General Public License v3 or later
# A copy of GNU AGPL v3 should have been included in this software package in LICENSE.txt.

""" NRPS/PKS analysis module
"""

from typing import List

from antismash.common import path
from antismash.config.args import ModuleArgs

from .html_output import will_handle, generate_sidepanel, generate_details_div, generate_js_domains
from .results import NRPS_PKS_Results
from .specific_analysis import specific_analysis

NAME = "nrps_pks"

SHORT_DESCRIPTION = "NRPS/PKS analysis"


def check_prereqs() -> List[str]:
    """ Check the prerequisites.
            java: NRPSPredictor, sandpuma
            muscle: sandpuma, at_analysis, kr_analysis, minowa, orderfinder
            hmmsearch: minowa
    """
    failure_messages = []
    for binary_name in ["hmmsearch", "muscle", "java"]:
        if path.locate_executable(binary_name) is None:
            failure_messages.append("Failed to locate executable for %r" %
                                    binary_name)

    return failure_messages


def get_arguments() -> ModuleArgs:
    """ Construct module arguments and for sub-modules """
    args = ModuleArgs('NRPS/PKS options', 'np', enabled_by_default=True)
    return args


def check_options(_options) -> List[str]:
    """ No options to check at the moment """
    return []


def regenerate_previous_results(json, record, _options) -> NRPS_PKS_Results:
    """ Regenerate results from previous run """
    return NRPS_PKS_Results.from_json(json, record)


def is_enabled(options):
    """ Whether the module is enabled """
    return not options.minimal or options.nrps_pks_enabled


def run_on_record(record, results, options):
    """ Analyse a record """
    if not results:
        results = NRPS_PKS_Results(record.id)
        specific_analysis(record, results, options)
    return results
