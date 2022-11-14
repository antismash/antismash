# License: GNU Affero General Public License v3 or later
# A copy of GNU AGPL v3 should have been included in this software package in LICENSE.txt.

""" Transcription factor binding site predictor.
"""

from typing import Any, Dict, List, Optional

from antismash.common.secmet import Record
from antismash.config import ConfigType
from antismash.config.args import ModuleArgs

from antismash.common import path

from .html_output import generate_html, generate_javascript_data, will_handle
from .tfbs_finder import PWM_PATH, run_tfbs_finder, TFBSFinderResults

NAME = "tfbs_finder"
SHORT_DESCRIPTION = "Detects transcription factor binding sites (TFBSs)"


def get_arguments() -> ModuleArgs:
    """ Builds any commandline argument constructs that may be required

        Returns:
            an empty or populated ModuleArgs instance
    """

    args = ModuleArgs('Transcription Factor Binding Site options', 'tfbs')

    args.add_analysis_toggle('--tfbs',
                             dest='tfbs',
                             default=False,
                             action='store_true',
                             help="Run TFBS finder on all gene clusters.")

    args.add_option('pvalue',
                    dest='pvalue',
                    type=float,
                    default=0.00001,
                    help="P-value for TFBS threshold setting (default: %(default)s).")

    args.add_option('range',
                    dest='range',
                    type=int,
                    default=50,
                    help=("The allowable overlap with gene start positions for "
                          "TFBSs in coding regions (default: %(default)s)."))
    return args


def check_options(options: ConfigType) -> List[str]:
    """ Checks that the provided options are compatible with each other

        Arguments:
            options: the current antismash config object

        Returns:
            a list of strings, each string being an issue with the given options
    """
    issues = []
    # MOODS doesn't output the p-value
    if options.tfbs_pvalue <= 0:
        issues.append(f"TFBS finder p-value is negative: {options.tfbs_pvalue}")
    if options.tfbs_range <= 0:
        issues.append(f"TFBS finder range is negative: {options.tfbs_range}")
    return issues


def check_prereqs(_options: ConfigType) -> List[str]:
    """ Check if all required applications are around """
    failure_messages = []

    if path.locate_file(PWM_PATH) is None:
        failure_messages.append("Failed to locate PWM file")

    return failure_messages


def is_enabled(options: ConfigType) -> bool:
    """ Returns True if the module is enabled with the options provided
    """
    return options.tfbs


def regenerate_previous_results(previous: Dict[str, Any], record: Record,
                                options: ConfigType) -> Optional[TFBSFinderResults]:
    """ Regenerate the previous results from JSON format.

        Arguments:
            previous: the previous results as from JSON
            record: the Record these previous results were originally created from
            options: the current antismash config object

        Returns:
            an instance of the module's ModuleResults implementation,
            or None if the current options require the analysis to be rerun or cannot be regenerated
    """
    if not previous:
        return None
    results = TFBSFinderResults.from_json(previous, record)
    if not results:
        return None
    if options.tfbs_pvalue < results.pvalue:
        return None
    if options.tfbs_range != results.start_overlap:
        return None
    return TFBSFinderResults.from_json(previous, record)


def run_on_record(record: Record, results: TFBSFinderResults, options: ConfigType
                  ) -> TFBSFinderResults:
    """ Run the analysis, unless the previous results apply to the given record

        Arguments:
            record: the Record being analysed
            results: an existing instance of the module's ModuleResults implementation (or None)
            options: the current antismash config object

        Returns:
            an instance of the module's ModuleResults implementation
    """
    if isinstance(results, TFBSFinderResults) and results.record_id == record.id:
        return results
    return run_tfbs_finder(record, options.tfbs_pvalue, options.tfbs_range)
