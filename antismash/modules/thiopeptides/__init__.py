# License: GNU Affero General Public License v3 or later
# A copy of GNU AGPL v3 should have been included in this software package in LICENSE.txt.

""" Thiopeptides detection module

"""

import os
import logging
from typing import Any, Dict, List, Optional

from antismash.common import comparippson, hmmer, path
from antismash.common.external.rodeo_svm.rebuild import pickle_classifier
from antismash.common.secmet import Record
from antismash.config import ConfigType
from antismash.config.args import ModuleArgs

from .specific_analysis import specific_analysis, ThioResults
from .html_output import generate_html, will_handle

NAME = "thiopeptides"
SHORT_DESCRIPTION = NAME.capitalize()


def check_options(_options: ConfigType) -> List[str]:
    """ No options to check for conflicts """
    return []


def is_enabled(options: ConfigType) -> bool:
    """ Only disabled if minimal and not re-enabled """
    return options.thiopeptides_enabled or not options.minimal


def get_arguments() -> ModuleArgs:
    """ Build arguments, no explicit options as enabled by default """
    args = ModuleArgs('Advanced options', 'thio', enabled_by_default=True)
    return args


def prepare_data(logging_only: bool = False) -> List[str]:
    """ Ensures packaged data is fully prepared

        Arguments:
            logging_only: whether to return error messages instead of raising exceptions

        Returns:
            a list of error messages (only if logging_only is True)
    """
    errors = comparippson.prepare_data(logging_only=logging_only)
    training_set = path.get_full_path(__file__, "data", "training_set.csv")
    expected = ["thiopeptide.scaler.pkl", "thiopeptide.classifier.pkl"]
    precursor_model = path.get_full_path(__file__, "data", "thiopep3.hmm")
    errors.extend(hmmer.ensure_database_pressed(precursor_model, return_not_raise=True))

    if all(os.path.exists(path.get_full_path(__file__, "data", filename)) for filename in expected):
        return errors
    try:
        pickle_classifier(training_set, prefix="thiopeptide", kernel='rbf', C=2.83e5, gamma=1e-9,
                          overwrite=not logging_only)
    except ValueError:
        if logging_only:
            errors.append("failed to rebuild thiopeptide classifier")
        else:
            raise
    return errors


def check_prereqs(options: ConfigType) -> List[str]:
    """ Check prereqs
            hmmpfam2: used to find extra HMM hits not in hmm_detection
    """
    failure_messages = []
    for binary_name in ['hmmpfam2']:
        if binary_name not in options.executables:
            failure_messages.append(f"Failed to locate executable for {binary_name}")
    failure_messages.extend(prepare_data(logging_only=True))
    failure_messages.extend(comparippson.check_prereqs(options))

    return failure_messages


def regenerate_previous_results(results: Dict[str, Any], record: Record, _options: ConfigType) -> Optional[ThioResults]:
    """ Regenerate previous results from JSON format """
    if not results:
        return None
    regenned = ThioResults.from_json(results, record)
    if not regenned:
        return None
    logging.debug("Reusing Thiopeptide results: %d clusters contained %d total motifs",
                  len(regenned.clusters_with_motifs), len(regenned.motifs))
    return regenned


def run_on_record(record: Record, previous_results: ThioResults, _options: ConfigType) -> ThioResults:
    """ Run the thiopeptide analysis on the given record"""
    if previous_results:
        return previous_results
    return specific_analysis(record)
