# License: GNU Affero General Public License v3 or later
# A copy of GNU AGPL v3 should have been included in this software package in LICENSE.txt.

""" The rule-based cluster detection section of clusterfinder """

import logging
from typing import Any, Dict, List, Optional

from antismash.common import path
from antismash.common.hmm_rule_parser import cluster_prediction
from antismash.common.module_results import DetectionResults
from antismash.common.secmet import Cluster, Record
from antismash.config import ConfigType
from antismash.config.args import ModuleArgs

NAME = "clusterfinder-rule-based"
SHORT_DESCRIPTION = "Less specific rule-based cluster detection"


def get_arguments() -> ModuleArgs:
    """ Sets up arguments for this module.
    """
    args = ModuleArgs('ClusterFinder rule-based options', 'cfrules')
    args.add_analysis_toggle('cfrules',
                             dest='cfrules',
                             action='store_true',
                             default=False,
                             help="Find loose rule-based clusters defined by cluster-finder.")
    return args


def is_enabled(options: ConfigType) -> bool:
    """  Uses the supplied options to determine if the module should be run
    """
    return options.cfrules


def check_options(_options: ConfigType) -> List[str]:
    """ Checks that extra options are valid """
    return []


class ClusterFinderRuleResults(DetectionResults):
    """ A container for clusters predicted by rules in this module """
    schema_version = 1

    def __init__(self, record_id: str, rule_results: cluster_prediction.RuleDetectionResults) -> None:
        super().__init__(record_id)
        self.rule_results = rule_results

    def to_json(self) -> Dict[str, Any]:
        return {"record_id": self.record_id,
                "schema_version": self.schema_version,
                "rule_results": self.rule_results.to_json()}

    @staticmethod
    def from_json(json: Dict[str, Any], record: Record) -> "ClusterFinderRuleResults":
        if json["schema_version"] != ClusterFinderRuleResults.schema_version:
            raise ValueError("Detection results have changed. No results can be reused.")
        assert json["record_id"] == record.id

        rule_parser_results = cluster_prediction.RuleDetectionResults.from_json(json["rule_results"], record)

        return ClusterFinderRuleResults(json["record_id"], rule_parser_results)

    def get_predicted_clusters(self) -> List[Cluster]:
        return self.rule_results.clusters


def check_prereqs() -> List[str]:
    """Don't check for prerequisites, we don't have any"""
    return []


def regenerate_previous_results(results: Dict[str, Any], record: Record,
                                _options: ConfigType) -> Optional[ClusterFinderRuleResults]:
    """ Regenerate previous results. """
    if not results:
        return None

    regenerated = ClusterFinderRuleResults.from_json(results, record)
    if not regenerated:
        return None
    for cluster in regenerated.get_predicted_clusters():
        record.add_cluster(cluster)
    return regenerated


def run_on_record(record: Record, results: Optional[ClusterFinderRuleResults],
                  _options: ConfigType) -> ClusterFinderRuleResults:
    """Load the data and run the cluster_predict tool"""
    if results:
        return results

    return ClusterFinderRuleResults(record.id, find_rule_based_clusters(record))


def find_rule_based_clusters(record: Record) -> cluster_prediction.RuleDetectionResults:
    """ Generate cluster predictions using HMM profile base rules.

        Arguments:
            record: the record to find clusters in

        Returns:
            RuleDetectionResults as provided by common.hmm_rule_parser.cluster_prediction
    """
    # TODO: don't go directly to hmm_detection's data
    signatures = path.get_full_path(__file__, "..", "hmm_detection", "data", "hmmdetails.txt")
    seeds = path.get_full_path(__file__, "..", "hmm_detection", "data", "bgc_seeds.hmm")
    rules = path.get_full_path(__file__, "cluster_rules.txt")
    equivalences = path.get_full_path(__file__, "..", "hmm_detection", "filterhmmdetails.txt")
    results = cluster_prediction.detect_clusters_and_signatures(record, signatures, seeds,
                                                                rules, equivalences, "cluster-finder")
    assert results is not None
    results.annotate_cds_features()
    logging.debug("ClusterFinder detected %d rule-based clusters", len(results.clusters))
    return results
