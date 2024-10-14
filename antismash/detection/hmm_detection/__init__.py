# License: GNU Affero General Public License v3 or later
# A copy of GNU AGPL v3 should have been included in this software package in LICENSE.txt.

""" General HMM detection module for detecting specific domains and defining
    clusters based on domains detected
"""

import importlib
import logging
import os
import pkgutil
from typing import Any, Dict, List, Iterable, Optional

from antismash.common import hmmer, path
from antismash.common.hmm_rule_parser import rule_parser
from antismash.common.hmm_rule_parser.cluster_prediction import (
    detect_protoclusters_and_signatures,
    RuleDetectionResults,
    Ruleset,
)
from antismash.common.hmm_rule_parser.structures import DynamicProfile, Multipliers
from antismash.common.module_results import DetectionResults
from antismash.common.secmet.record import Record
from antismash.common.secmet.features import Protocluster
from antismash.config import ConfigType
from antismash.config.args import ModuleArgs, SplitCommaAction
from antismash.detection import DetectionStage
from antismash.detection.hmm_detection.signatures import get_signature_profiles
from antismash.detection.hmm_detection.categories import get_rule_categories

NAME = "hmmdetection"
SHORT_DESCRIPTION = "HMM signature detection"
DETECTION_STAGE = DetectionStage.AREA_FORMATION

SIGNATURE_FILE = path.get_full_path(__file__, "data", "hmmdetails.txt")
HMM_FILE = path.get_full_path(__file__, "data", "bgc_seeds.hmm")
CATEGORIES = {cat.name for cat in get_rule_categories()}
EQUIVALENCE_GROUPS = path.get_full_path(__file__, "filterhmmdetails.txt")

_STRICTNESS_LEVELS = ["strict", "relaxed", "loose"]

_RULESETS: dict[tuple[str, tuple[str, ...], tuple[str, ...], Multipliers], Ruleset] = {}


def _get_dynamic_profiles() -> Dict[str, DynamicProfile]:
    """ Gather all the dynamic profiles """
    profiles = {}
    for module_data in pkgutil.walk_packages([path.get_full_path(__file__, "dynamic_profiles")]):
        module = importlib.import_module(f"antismash.detection.hmm_detection.dynamic_profiles.{module_data.name}")
        contains_profiles = False
        for name, profile in vars(module).items():
            if not isinstance(profile, DynamicProfile):
                continue
            if profile.name in profiles:
                raise ValueError(f"duplicate dynamic profile detected: dynamic_profiles.{module_data.name}.{name}")
            profiles[profile.name] = profile
            contains_profiles = True
        if not contains_profiles:
            raise ValueError(f"dynamic profile subpackage {module_data.name} has no DynamicProfile instances")
    return profiles


DYNAMIC_PROFILES = _get_dynamic_profiles()


def _get_rule_files_for_strictness(strictness: str) -> List[str]:
    """ Returns a list of appropriate rule files for the given strictness level """
    assert strictness in _STRICTNESS_LEVELS, strictness
    files = []
    for level in _STRICTNESS_LEVELS[:_STRICTNESS_LEVELS.index(strictness) + 1]:
        files.append(path.get_full_path(__file__, "cluster_rules", f"{level}.txt"))
    return files


def get_ruleset(options: ConfigType) -> Ruleset:
    """ Builds a Ruleset instance configured to match the provided options

        Arguments:
            options: the antiSMASH config object

        Returns:
            a Ruleset instance
    """
    strictness = options.hmmdetection_strictness
    name_subset = set(options.hmmdetection_limit_to_rules)
    category_subset = set(options.hmmdetection_limit_to_categories)
    multipliers = Multipliers()
    if options.taxon == "fungi":
        multipliers = Multipliers(
            options.hmmdetection_fungal_cutoff_multiplier,
            options.hmmdetection_fungal_neighbourhood_multiplier,
        )
    # the cache key needs to be immutable
    key = (strictness, tuple(name_subset), tuple(category_subset), multipliers)

    # return any existing ruleset
    ruleset = _RULESETS.get(key)
    if ruleset:
        return ruleset

    # otherwise make a default ruleset for the strictness
    ruleset = Ruleset.from_files(SIGNATURE_FILE, HMM_FILE, _get_rule_files_for_strictness(strictness),
                                 CATEGORIES, EQUIVALENCE_GROUPS, "rule-based-clusters",
                                 dynamic_profiles=DYNAMIC_PROFILES)

    # limit the rules used, if relevant
    rules: Iterable[rule_parser.DetectionRule] = ruleset.rules
    if name_subset:
        rules = filter(lambda rule: rule.name in name_subset, rules)
    if category_subset:
        rules = filter(lambda rule: rule.category in category_subset, rules)

    ruleset = ruleset.copy_with_replacements(rules=list(rules), multipliers=multipliers)

    # update the cache
    _RULESETS[key] = ruleset

    return ruleset


class HMMDetectionResults(DetectionResults):
    """ A container for clusters predicted by rules in this module """
    schema_version = 2

    def __init__(self, record_id: str, rule_results: RuleDetectionResults, enabled_types: List[str],
                 strictness: str) -> None:
        super().__init__(record_id)
        self.rule_results = rule_results
        self.enabled_types = enabled_types
        if strictness not in _STRICTNESS_LEVELS:
            raise ValueError(f"unknown strictness level: {strictness}")
        self.strictness = strictness

    def to_json(self) -> Dict[str, Any]:
        return {"record_id": self.record_id,
                "schema_version": self.schema_version,
                "enabled_types": self.enabled_types,
                "rule_results": self.rule_results.to_json(),
                "strictness": self.strictness}

    @staticmethod
    def from_json(json: Dict[str, Any], record: Record) -> "HMMDetectionResults":
        if json["schema_version"] != HMMDetectionResults.schema_version:
            raise ValueError("Detection results have changed. No results can be reused.")
        assert json["record_id"] == record.id

        rule_results = RuleDetectionResults.from_json(json["rule_results"], record)
        if rule_results is None:
            raise ValueError("Detection results have changed. No results can be reused")

        return HMMDetectionResults(json["record_id"], rule_results, json["enabled_types"],
                                   json.get("strictness", "relaxed"))

    def get_predicted_protoclusters(self) -> List[Protocluster]:
        return self.rule_results.protoclusters


def _get_rules(strictness: str) -> list[rule_parser.DetectionRule]:
    signature_names = {sig.name for sig in get_signature_profiles()}
    signature_names.update(set(DYNAMIC_PROFILES))
    category_names = {cat.name for cat in get_rule_categories()}
    rules: List[rule_parser.DetectionRule] = []
    aliases: Dict[str, List[rule_parser.Token]] = {}
    for rule_file in _get_rule_files_for_strictness(strictness):
        with open(rule_file, encoding="utf-8") as rulefile:
            rules = rule_parser.Parser("".join(rulefile.readlines()), signature_names,
                                       category_names, rules, aliases).rules
    return rules


def get_arguments() -> ModuleArgs:
    """ Constructs commandline arguments and options for this module
    """
    args = ModuleArgs('HMM detection options', 'hmmdetection')
    args.add_option('strictness',
                    dest='strictness',
                    type=str,
                    choices=["strict", "relaxed", "loose"],
                    default="relaxed",
                    help=("Defines which level of strictness to use for "
                          "HMM-based cluster detection, (default: %(default)s)."))
    args.add_option("fungal-cutoff-multiplier",
                    dest="fungal_cutoff_multiplier",
                    type=float,
                    default=1.0,
                    help=("Sets the multiplier for rule cutoffs in fungal inputs "
                          "(default: %(default)s)."))
    args.add_option("fungal-neighbourhood-multiplier",
                    dest="fungal_neighbourhood_multiplier",
                    type=float,
                    default=1.5,
                    help=("Sets the multiplier for rule neighbourhoods in fungal "
                          "inputs (default: %(default)s)."))
    args.add_option("limit-to-rule-names",
                    dest="limit_to_rules",
                    metavar="RULE1[,RULE2,...]",
                    action=SplitCommaAction,
                    default=[],
                    help="Restrict detection to the named rules (default: no limits).")
    args.add_option("limit-to-rule-categories",
                    dest="limit_to_categories",
                    metavar="CATEGORY1[,CATEGORY2,...]",
                    action=SplitCommaAction,
                    default=[],
                    help="Restrict detection to the given rules (default: no limits).")

    return args


def check_options(options: ConfigType) -> List[str]:
    """ Checks the options to see if there are any issues before
        running any analyses
    """
    issues = []
    if options.hmmdetection_strictness not in _STRICTNESS_LEVELS:
        issues.append(f"Unknown strictness level: {options.strictness}")
    cutoff = options.hmmdetection_fungal_cutoff_multiplier
    if cutoff <= 0:
        issues.append(f"Invalid fungal cutoff multiplier: {cutoff}")
    neighbourhood = options.hmmdetection_fungal_neighbourhood_multiplier
    if neighbourhood <= 0:
        issues.append(f"Invalid fungal neighbourhood multiplier: {neighbourhood}")

    all_rule_names = set(rule.name for rule in _get_rules(options.hmmdetection_strictness))

    name_subset = set(options.hmmdetection_limit_to_rules)
    if name_subset:
        unknown = name_subset - all_rule_names
        if unknown:
            issues.append(f"Unknown rules in requested rule subset: {unknown}")
    category_subset = set(options.hmmdetection_limit_to_categories)
    if category_subset:
        unknown = category_subset - CATEGORIES
        if unknown:
            issues.append(f"Unknown rules in requested rule category subset: {unknown}")

    # invalid multipliers will cause ruleset creation to fail for the same reasons
    # so don't continue into testing ruleset creation
    if issues:
        return issues

    try:
        get_ruleset(options)
    except ValueError as err:
        issues.append(str(err))

    return issues


def is_enabled(_options: ConfigType) -> bool:
    """  Uses the supplied options to determine if the module should be run
    """
    # in this case, yes, always
    return True


def regenerate_previous_results(results: Dict[str, Any], record: Record,
                                options: ConfigType) -> Optional[HMMDetectionResults]:
    """ Regenerate previous results. """
    if not results:
        return None
    regenerated = HMMDetectionResults.from_json(results, record)
    if regenerated.strictness != options.hmmdetection_strictness:
        logging.warning("Ignoring hmmdetection strictness option %r, reusing %r from results",
                        options.hmmdetection_strictness, regenerated.strictness)
    if set(regenerated.enabled_types) != get_ruleset(options).get_rule_names():
        raise RuntimeError("Protocluster types supported by HMM detection have changed, all results invalid")
    if options.taxon == "fungi":
        if regenerated.rule_results.multipliers.cutoff != options.hmmdetection_fungal_cutoff_multiplier:
            raise RuntimeError("Protocluster cutoff multiplier changed, previous results are incompatible")
        if regenerated.rule_results.multipliers.neighbourhood != options.hmmdetection_fungal_neighbourhood_multiplier:
            raise RuntimeError("Protocluster neighbourhood multiplier changed, previous results are incompatible")
    regenerated.rule_results.annotate_cds_features()
    return regenerated


def run_on_record(record: Record, previous_results: Optional[HMMDetectionResults],
                  options: ConfigType) -> HMMDetectionResults:
    """ Runs hmm_detection on the provided record.
    """
    if previous_results:
        return previous_results

    strictness = options.hmmdetection_strictness
    logging.info("HMM detection using strictness: %s", strictness)

    ruleset = get_ruleset(options)
    cluster_types = list(ruleset.get_rule_names())
    if options.hmmdetection_limit_to_rules or options.hmmdetection_limit_to_categories:
        logging.info("HMM detection restricted to these rules: %s", cluster_types)

    results = detect_protoclusters_and_signatures(record, ruleset)
    results.annotate_cds_features()
    return HMMDetectionResults(record.id, results, cluster_types, strictness)


def prepare_data(logging_only: bool = False) -> List[str]:
    """ Ensures packaged data is fully prepared

        Arguments:
            logging_only: whether to return error messages instead of raising exceptions

        Returns:
            a list of error messages (only if logging_only is True)
    """
    failure_messages = []

    # Check that hmmdetails.txt is readable and well-formatted
    try:
        profiles = get_signature_profiles()
    except ValueError as err:
        if not logging_only:
            raise
        return [str(err)]

    # the path to the markov model
    seeds_hmm = path.get_full_path(__file__, 'data', 'bgc_seeds.hmm')
    hmm_files = [os.path.join("data", sig.hmm_file) for sig in profiles]
    # include the listing, since tools like wget will keep modified timestamps on the HMMs
    description_file = path.get_full_path(__file__, 'data', 'hmmdetails.txt')
    outdated = False
    if not path.locate_file(seeds_hmm):
        logging.debug("%s: %s doesn't exist, regenerating", NAME, seeds_hmm)
        outdated = True
    else:
        seeds_timestamp = os.path.getmtime(seeds_hmm)
        for component in hmm_files + [description_file]:
            if os.path.getmtime(component) > seeds_timestamp:
                logging.debug("%s out of date, regenerating", seeds_hmm)
                outdated = True
                break

    # regenerate if missing or out of date
    if outdated:
        # try to generate file from all specified profiles in hmmdetails
        try:
            with open(seeds_hmm, "w", encoding="utf-8") as all_hmms_handle:
                for hmm_file in hmm_files:
                    with open(path.get_full_path(__file__, hmm_file), "r", encoding="utf-8") as handle:
                        all_hmms_handle.write(handle.read())
        except OSError:
            if not logging_only:
                raise
            failure_messages.append(f"Failed to generate file {seeds_hmm!r}")

    # if regeneration failed, don't try to run hmmpress
    if failure_messages:
        return failure_messages

    failure_messages.extend(hmmer.ensure_database_pressed(seeds_hmm, return_not_raise=logging_only))

    return failure_messages


def check_prereqs(options: ConfigType) -> List[str]:
    """ Check that prereqs are satisfied. hmmpress is only required if the
        databases have not yet been generated.
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
