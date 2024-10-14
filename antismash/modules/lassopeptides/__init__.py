# License: GNU Affero General Public License v3 or later
# A copy of GNU AGPL v3 should have been included in this software package in LICENSE.txt.

"""Lassopeptides detection module

"""

import os
from typing import Any, Dict, List, Optional

from antismash.common import comparippson, path
from antismash.common.external.rodeo_svm.rebuild import pickle_classifier
from antismash.common.secmet import Record
from antismash.config import ConfigType
from antismash.config.args import ModuleArgs

from .config import get_config as local_config
from .specific_analysis import specific_analysis as run_analysis, LassoResults
from .html_output import generate_html, will_handle

NAME = "lassopeptides"
SHORT_DESCRIPTION = "lassopeptide precursor prediction"


def prepare_data(logging_only: bool = False) -> List[str]:
    """ Ensures packaged data is fully prepared

        Arguments:
            logging_only: whether to return error messages instead of raising exceptions

        Returns:
            a list of error messages (only if logging_only is True)
    """
    errors = comparippson.prepare_data(logging_only=logging_only)
    training_set = path.get_full_path(__file__, "data", "training_set.csv")
    expected = ["lassopeptide.scaler.pkl", "lassopeptide.classifier.pkl"]
    if all(os.path.exists(path.get_full_path(__file__, "data", filename)) for filename in expected):
        return errors
    try:
        pickle_classifier(training_set, prefix="lassopeptide", kernel='rbf', C=2.83e5, gamma=1e-8,
                          overwrite=not logging_only)
    except ValueError:
        if logging_only:
            errors.append("failed to rebuild lassopeptide classifier")
        else:
            raise
    return errors


def check_prereqs(options: ConfigType) -> List[str]:
    """ Checks if the required external programs are available """
    failure_messages = []
    for binary_name, optional in [('hmmpfam2', False), ('fimo', True)]:
        present = True
        if binary_name not in options.executables:
            present = False
            if not optional:
                failure_messages.append(f"Failed to locate executable for {binary_name}")
        if binary_name == "fimo":
            local_config().fimo_present = present

    failure_messages.extend(prepare_data(logging_only=True))
    failure_messages.extend(comparippson.check_prereqs(options))
    return failure_messages


def get_arguments() -> ModuleArgs:
    """ Runs by default, but add minimal's --enable option """
    args = ModuleArgs('Advanced options', 'lasso', enabled_by_default=True)
    return args


def check_options(_options: ConfigType) -> List[str]:
    """ No options here to check, so just return """
    return []


def is_enabled(options: ConfigType) -> bool:
    """ Will the module run with the given options """
    return not options.minimal or options.lassopeptides_enabled


def regenerate_previous_results(previous: Dict[str, Any], record: Record,
                                _options: ConfigType) -> Optional[LassoResults]:
    """ Regenerate a results object from the given data """
    return LassoResults.from_json(previous, record)


def run_on_record(record: Record, results: LassoResults, _options: ConfigType) -> LassoResults:
    """ Finds all precursors within lassopeptide clusters """
    if results and isinstance(results, LassoResults):
        return results
    return run_analysis(record)
