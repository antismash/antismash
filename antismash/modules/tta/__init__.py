# License: GNU Affero General Public License v3 or later
# A copy of GNU AGPL v3 should have been included in this software package in LICENSE.txt.

"""Identify TTA codons in BGCs"""

from typing import Any, Dict, List, Optional

from antismash.common.secmet import Record
from antismash.config import ConfigType
from antismash.config.args import ModuleArgs
from antismash.modules.tta.tta import detect, TTAResults

NAME = "tta"
SHORT_DESCRIPTION = "TTA detection"


def get_arguments() -> ModuleArgs:
    """ Build and return arguments. No extra options beyond a switch to enable """
    args = ModuleArgs('Additional analysis', 'tta')
    args.add_analysis_toggle('--tta',
                             dest='tta',
                             action='store_true',
                             default=False,
                             help="Run TTA codon detection module.")
    return args


def check_options(_options: ConfigType) -> List[str]:
    """ Checks options for conflicts.
        No extra options, so they can't have conflicts.
    """
    return []


def check_prereqs() -> List[str]:
    """Check for prerequisites"""
    # No external dependencies
    return []


def is_enabled(options: ConfigType) -> bool:
    """ Should the module be run with these options """
    return options.tta


def regenerate_previous_results(previous: Dict[str, Any], record: Record,
                                _options: ConfigType) -> Optional[TTAResults]:
    """ Regenerate the previous results from JSON format. """
    if not previous:
        return None
    return TTAResults.from_json(previous, record)


def run_on_record(record: Record, results: TTAResults, options: ConfigType) -> TTAResults:
    """ Run the analysis, unless the previous results apply to the given record """
    if isinstance(results, TTAResults) and results.record_id == record.id:
        return results
    return detect(record, options)
