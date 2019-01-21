# License: GNU Affero General Public License v3 or later
# A copy of GNU AGPL v3 should have been included in this software package in LICENSE.txt.

""" NRPS/PKS analysis module
"""

from typing import Any, Dict, List, Optional

from antismash.common import path
from antismash.common.secmet import Record
from antismash.config import ConfigType
from antismash.config.args import ModuleArgs

from .html_output import generate_html, will_handle
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


def check_options(_options: ConfigType) -> List[str]:
    """ No options to check at the moment """
    return []


def regenerate_previous_results(json: Dict[str, Any], record: Record,
                                _options: ConfigType) -> Optional[NRPS_PKS_Results]:
    """ Regenerate results from previous run """
    return NRPS_PKS_Results.from_json(json, record)


def is_enabled(options: ConfigType) -> bool:
    """ Whether the module is enabled """
    return not options.minimal or options.nrps_pks_enabled


def run_on_record(record: Record, results: Optional[NRPS_PKS_Results], options: ConfigType) -> NRPS_PKS_Results:
    """ Analyse a record """
    if not results:
        results = NRPS_PKS_Results(record.id)
        specific_analysis(record, results, options)
    return results
