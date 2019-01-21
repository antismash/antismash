# License: GNU Affero General Public License v3 or later
# A copy of GNU AGPL v3 should have been included in this software package in LICENSE.txt.

"""Sactipeptides detection module

"""

from typing import Any, Dict, List, Optional

from antismash.common import path
from antismash.common.secmet import Record
from antismash.config import ConfigType
from antismash.config.args import ModuleArgs

from .specific_analysis import specific_analysis, SactiResults
from .html_output import generate_html, will_handle

NAME = "sactipeptides"
SHORT_DESCRIPTION = "sactipeptide detection"


def check_prereqs() -> List[str]:
    """ Ensures all required external programs are available """
    failure_messages = []
    for binary_name in ['hmmpfam2']:
        if not path.locate_executable(binary_name):
            failure_messages.append("Failed to locate executable for %r" %
                                    binary_name)
    return failure_messages


def get_arguments() -> ModuleArgs:
    """ Runs by default, but add minimal's --enable option """
    args = ModuleArgs('Advanced options', 'sacti', enabled_by_default=True)
    return args


def check_options(_options: ConfigType) -> List[str]:
    """ No specific options to check """
    return []


def is_enabled(options: ConfigType) -> bool:
    """ Should the module run with the given options """
    return not options.minimal or options.sactipeptides_enabled


def regenerate_previous_results(previous: Dict[Any, str], record: Record,
                                _options: ConfigType) -> Optional[SactiResults]:
    """ Regenerate a SactiResults object from a JSON-like dict """
    return SactiResults.from_json(previous, record)


def run_on_record(record: Record, results: Optional[SactiResults], _options: ConfigType) -> SactiResults:
    """ Run the sactipeptide analysis over the given record """
    if results:
        return results
    return specific_analysis(record)
