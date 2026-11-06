# License: GNU Affero General Public License v3 or later
# A copy of GNU AGPL v3 should have been included in this software package in LICENSE.txt.

""" Detection of subclusters, i.e. sets of genes responsible for producing a
    specific chemical moiety
"""
import logging
import os
from typing import Any, Optional

from antismash.common import hmmer, path
from antismash.common.hmm_rule_parser.cluster_prediction import detect_protoclusters_and_signatures
from antismash.common.secmet import Record
from antismash.config import ConfigType
from antismash.config.args import ModuleArgs
from antismash.detection import DetectionStage

from .predictions import build_predictions
from .results import SubclusterDetectionResults
from .ruleset import STRICTNESS_LEVELS, get_ruleset
from .signatures import AGGREGATE_HMM_FILE, DETAILS_FILE, get_signatures
from .subregions import SubRegionMode, build_subregions
from .html_output import generate_html, will_handle, generate_javascript_data  # noqa: F401

NAME = "subclusters"
SHORT_DESCRIPTION = "Detects subclusters for specific chemical moieties, e.g. precursor molecules."
DETECTION_STAGE = DetectionStage.AREA_FORMATION


def get_arguments() -> ModuleArgs:
    """ Constructs commandline arguments and options for this module

        Returns:
            an empty or populated set of module arguments
    """
    args = ModuleArgs("Subcluster detection options", "subclusters")

    args.add_analysis_toggle('--subclusters',
                             dest='subclusters',
                             default=False,
                             action='store_true',
                             help="Run Subcluster detection.")
    args.add_option('strictness',
                    dest='strictness',
                    type=str,
                    choices=list(STRICTNESS_LEVELS),
                    default="relaxed",
                    help=("Defines which level of strictness to use for "
                          "subcluster detection. Levels are cumulative, so "
                          "looser levels also include all stricter rules "
                          "(default: %(default)s)."))
    args.add_option('region-mode',
                    dest='subregion_mode',
                    type=str,
                    choices=list(SubRegionMode),
                    default=SubRegionMode.CREATE,
                    help=("How far detected subclusters may influence region formation. "
                          "The 'clip' option keeps only the part of a subcluster "
                          "overlapping an area found by another detection module, "
                          "truncated to that area, so regions can never grow. The "
                          "'extend' option keeps overlapping subclusters in full, so "
                          "they can extend an existing region. The 'create' option keeps "
                          "every subcluster in full, so they can also create "
                          "standalone regions. (default: %(default)s)."))
    return args


def check_options(options: ConfigType) -> list[str]:
    """ Checks the options to see if there are any issues before
        running any analyses
    """
    errors = []
    if options.subclusters_strictness not in STRICTNESS_LEVELS:
        errors.append(f"Unknown subcluster strictness level: {options.subclusters_strictness}")

    try:
        SubRegionMode(options.subclusters_subregion_mode)
    except ValueError:
        errors.append(f"Unknown subcluster region mode: {options.subclusters_subregion_mode}")

    return errors


def is_enabled(options: ConfigType) -> bool:
    """ Uses the supplied options to determine if the module should be run """
    return options.subclusters


def prepare_data(logging_only: bool = False) -> list[str]:
    """ Ensures packaged data is fully prepared

        Aggregates the individual subcluster HMM signatures into a single
        combined profile database and presses it with hmmpress, regenerating
        whenever it is missing or out of date.

        Arguments:
            logging_only: whether to return error messages instead of raising exceptions

        Returns:
            a list of error messages (only if logging_only is True)
    """
    failure_messages: list[str] = []

    # check that the signature details are readable and well-formatted
    try:
        signatures = get_signatures()
    except ValueError as err:
        if not logging_only:
            raise
        return [str(err)]
    hmm_files = [signature.hmm_file for signature in signatures.values()]

    # aggregate_profiles only compares the aggregate against the profiles it was built
    # from, so changes to the details describing them have to be caught here instead
    force_replace = not (path.locate_file(AGGREGATE_HMM_FILE)
                         and os.path.getmtime(DETAILS_FILE) < os.path.getmtime(AGGREGATE_HMM_FILE))

    failure_messages.extend(hmmer.aggregate_profiles(
        AGGREGATE_HMM_FILE, hmm_files, force_replace=force_replace, return_not_raise=logging_only)
        )

    return failure_messages


def check_prereqs(options: ConfigType) -> list[str]:
    """ Checks that the required executables are available and that the
        packaged data is prepared
    """
    failure_messages = []
    for binary_name in ["hmmsearch", "hmmpress"]:
        if binary_name not in options.executables:
            failure_messages.append(f"Failed to locate executable for {binary_name!r}")

    # no point checking the data if we can't use it
    if failure_messages:
        return failure_messages

    failure_messages.extend(prepare_data(logging_only=True))

    return failure_messages


def regenerate_previous_results(results: dict[str, Any], record: Record,
                                options: ConfigType) -> Optional[SubclusterDetectionResults]:
    """ Regenerates previous results, if they are still compatible with the
        current options

        Since regions are rebuilt from scratch when reusing results, any change
        that could alter region formation would leave other modules' results
        incompatible.

        Arguments:
            results: the JSON representation of the previous results
            record: the record the results were generated for
            options: the current antiSMASH config

        Returns:
            the regenerated results, or None if there are no previous results
            or they were discarded so that detection is rerun

        Raises:
            RuntimeError: if the previous results are incompatible and rerunning
                          detection could alter regions
    """
    if not results:
        return None

    current_mode = SubRegionMode(options.subclusters_subregion_mode)
    current_rule_names = get_ruleset(options.subclusters_strictness).get_rule_names()

    mode_changed = results.get("subregion_mode") != current_mode
    rules_changed = set(results.get("rule_names", [])) != current_rule_names

    # full antismash rerun necessary if the mode changes
    if mode_changed:
        raise RuntimeError("Subcluster region mode changed, no results can be reused")

    # rule set must be unchanged to regenerate results
    if rules_changed:
        regenerated = None
    else:
        # from_json returns None if the results format has changed
        regenerated = SubclusterDetectionResults.from_json(results, record)

    # only in CLIP mode it's safe to rerun subcluster module without rerunning everything
    if regenerated is None:
        if current_mode == SubRegionMode.CLIP:
            logging.debug("Subcluster rules or results format changed, "
                          "discarding previous subcluster results")
        else:
            raise RuntimeError(
                "Subcluster rules or results format changed, no results can be reused")

    return regenerated


def run_on_record(record: Record, previous_results: Optional[SubclusterDetectionResults],
                  options: ConfigType) -> SubclusterDetectionResults:
    """ Runs subcluster detection over the given record """
    if previous_results:
        return previous_results

    current_strictness = options.subclusters_strictness
    ruleset = get_ruleset(current_strictness)
    rule_results = detect_protoclusters_and_signatures(record, ruleset)
    # the shared detection labels protoclusters as "rule-based-clusters", so use this module's name
    for protocluster in rule_results.protoclusters:
        protocluster.tool = rule_results.tool

    predictions = build_predictions(rule_results, current_strictness)

    subregion_mode = SubRegionMode(options.subclusters_subregion_mode)
    subregions = build_subregions(record, rule_results.protoclusters,
                                  tool=rule_results.tool, mode=subregion_mode)

    return SubclusterDetectionResults(
        record_id=record.id,
        rule_results=rule_results,
        rule_names=ruleset.get_rule_names(),
        strictness=current_strictness,
        subregion_mode=subregion_mode,
        predictions=predictions,
        subregions=subregions,
    )
