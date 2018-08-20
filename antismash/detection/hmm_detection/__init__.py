# License: GNU Affero General Public License v3 or later
# A copy of GNU AGPL v3 should have been included in this software package in LICENSE.txt.

""" General HMM detection module for detecting specific domains and defining
    clusters based on domains detected
"""

import os
from typing import Any, Dict, List, Optional

import antismash.common.path as path
from antismash.common.subprocessing import run_hmmpress

from antismash.config import ConfigType
from antismash.config.args import ModuleArgs
from antismash.common.hmm_rule_parser import rule_parser
from antismash.common.hmm_rule_parser.cluster_prediction import detect_clusters_and_signatures, RuleDetectionResults
from antismash.common.module_results import DetectionResults
from antismash.common.secmet.record import Record
from antismash.common.secmet.features import Cluster
from antismash.detection.hmm_detection.signatures import get_signature_profiles

NAME = "hmmdetection"
SHORT_DESCRIPTION = "HMM signature detection"


class HMMDetectionResults(DetectionResults):
    """ A container for clusters predicted by rules in this module """
    schema_version = 1

    def __init__(self, record_id: str, rule_results: RuleDetectionResults, enabled_types: List[str]) -> None:
        super().__init__(record_id)
        self.rule_results = rule_results
        self.enabled_types = enabled_types

    def to_json(self) -> Dict[str, Any]:
        return {"record_id": self.record_id,
                "schema_version": self.schema_version,
                "enabled_types": self.enabled_types,
                "rule_results": self.rule_results.to_json()}

    @staticmethod
    def from_json(json: Dict[str, Any], record: Record) -> "HMMDetectionResults":
        if json["schema_version"] != HMMDetectionResults.schema_version:
            raise ValueError("Detection results have changed. No results can be reused.")
        assert json["record_id"] == record.id

        return HMMDetectionResults(json["record_id"], RuleDetectionResults.from_json(json["rule_results"], record),
                                   json["enabled_types"])

    def get_predicted_clusters(self) -> List[Cluster]:
        return self.rule_results.clusters


def get_supported_cluster_types() -> List[str]:
    """ Returns a list of all cluster types for which there are rules
    """
    signature_names = {sig.name for sig in get_signature_profiles()}
    with open(path.get_full_path(__file__, 'cluster_rules.txt'), "r") as rulefile:
        rules = rule_parser.Parser("".join(rulefile.readlines()), signature_names).rules
        clustertypes = [rule.name for rule in rules]
    return clustertypes


def get_arguments() -> ModuleArgs:
    """ Constructs commandline arguments and options for this module
    """
    return ModuleArgs('Advanced options', 'hmmdetection')


def check_options(_options: ConfigType) -> List[str]:
    """ Checks the options to see if there are any issues before
        running any analyses
    """
    return []


def is_enabled(_options: ConfigType) -> bool:
    """  Uses the supplied options to determine if the module should be run
    """
    # in this case, yes, always
    return True


def regenerate_previous_results(results: Dict[str, Any], record: Record,
                                _options: ConfigType) -> Optional[HMMDetectionResults]:
    """ Regenerate previous results. """
    if not results:
        return None
    regenerated = HMMDetectionResults.from_json(results, record)
    if set(regenerated.enabled_types) != set(get_supported_cluster_types()):
        raise RuntimeError("Cluster types supported by HMM detection have changed, all results invalid")
    regenerated.rule_results.annotate_cds_features()
    return regenerated


def run_on_record(record: Record, previous_results: Optional[HMMDetectionResults],
                  _options: ConfigType) -> HMMDetectionResults:
    """ Runs hmm_detection on the provided record.
    """
    if previous_results:
        return previous_results

    signatures = path.get_full_path(__file__, "data", "hmmdetails.txt")
    seeds = path.get_full_path(__file__, "data", "bgc_seeds.hmm")
    rules = path.get_full_path(__file__, "cluster_rules.txt")
    equivalences = path.get_full_path(__file__, "filterhmmdetails.txt")
    results = detect_clusters_and_signatures(record, signatures, seeds, rules, equivalences,
                                             "rule-based-clusters")
    results.annotate_cds_features()
    return HMMDetectionResults(record.id, results, get_supported_cluster_types())


def check_prereqs() -> List[str]:
    """ Check that prereqs are satisfied. hmmpress is only required if the
        databases have not yet been generated.
    """
    failure_messages = []
    for binary_name, optional in [('hmmsearch', False), ('hmmpress', False)]:
        if path.locate_executable(binary_name) is None and not optional:
            failure_messages.append("Failed to locate executable for %r" %
                                    binary_name)

    profiles = None
    # Check that hmmdetails.txt is readable and well-formatted
    try:
        profiles = get_signature_profiles()
    except ValueError as err:
        failure_messages.append(str(err))

    # the path to the markov model
    hmm = path.get_full_path(__file__, 'data', 'bgc_seeds.hmm')
    hmm_files = [os.path.join("data", sig.hmm_file) for sig in profiles]
    if path.locate_file(hmm) is None:
        # try to generate file from all specified profiles in hmmdetails
        try:
            with open(hmm, 'w') as all_hmms_handle:
                for hmm_file in hmm_files:
                    with open(path.get_full_path(__file__, hmm_file), 'r') as handle:
                        all_hmms_handle.write(handle.read())
        except OSError:
            failure_messages.append('Failed to generate file {!r}'.format(hmm))

    # if previous steps have failed, the remainder will too, so don't try
    if failure_messages:
        return failure_messages

    binary_extensions = ['.h3f', '.h3i', '.h3m', '.h3p']
    for ext in binary_extensions:
        binary = "{}{}".format(hmm, ext)
        if path.locate_file(binary) is None:
            result = run_hmmpress(hmm)
            if not result.successful():
                failure_messages.append('Failed to hmmpress {!r}: {}'.format(hmm, result.stderr))
            break

    return failure_messages
