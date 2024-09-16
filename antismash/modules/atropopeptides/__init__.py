# License: GNU Affero General Public License v3 or later
# A copy of GNU AGPL v3 should have been included in this software package in LICENSE.txt.

"""Atropopeptides detection module

"""

from typing import Any, Optional

from antismash.common import comparippson, path
from antismash.common.secmet import Record
from antismash.config import ConfigType
from antismash.config.args import ModuleArgs

from .specific_analysis import specific_analysis as run_analysis, AtropoResults
from .html_output import generate_html, will_handle


NAME = "atropopeptides"
SHORT_DESCRIPTION = "atropopeptide precursor prediction"


def prepare_data(logging_only: bool = False) -> list[str]:
    """ Ensures packaged data is fully prepared

        Arguments:
            logging_only: whether to return error messages instead of raising exceptions

        Returns:
            a list of error messages (only if logging_only is True)
    """
    errors = comparippson.prepare_data(logging_only=logging_only)
    return errors


def check_prereqs(options: ConfigType) -> list[str]:
    """ Checks if the required external programs are available """
    failure_messages = []
    failure_messages.extend(prepare_data(logging_only=True))
    failure_messages.extend(comparippson.check_prereqs(options))
    return failure_messages


def get_arguments() -> ModuleArgs:
    """ Runs by default, but add minimal's --enable option """
    args = ModuleArgs('Advanced options', 'lasso', enabled_by_default=True)
    return args


def check_options(_options: ConfigType) -> list[str]:
    """ No options here to check, so just return """
    return []


def is_enabled(options: ConfigType) -> bool:
    """ Will the module run with the given options """
    return not options.minimal or options.lassopeptides_enabled


def regenerate_previous_results(previous: dict[str, Any], record: Record,
                                _options: ConfigType) -> Optional[AtropoResults]:
    """ Regenerate a results object from the given data """
    return AtropoResults.from_json(previous, record)


def run_on_record(record: Record, results: AtropoResults, _options: ConfigType) -> AtropoResults:
    """ Finds all precursors within lassopeptide clusters """
    if results and isinstance(results, AtropoResults):
        return results
    return run_analysis(record)

