# License: GNU Affero General Public License v3 or later
# A copy of GNU AGPL v3 should have been included in this software package in LICENSE.txt.

""" Probability based cluster detection
"""

import logging
from typing import Any, Dict, List, Optional

from antismash.common import path
from antismash.common.hmm_rule_parser.cluster_prediction import detect_borders_and_signatures
from antismash.common.module_results import DetectionResults
from antismash.common.secmet import ClusterBorder, Record
from antismash.config import ConfigType
from antismash.config.args import ModuleArgs

from .probabilistic import find_probabilistic_clusters, get_pfam_probabilities

NAME = "clusterfinder"
SHORT_DESCRIPTION = "Statistics-based cluster detection"

PUTATIVE_PRODUCT = "cf_putative"


def get_arguments() -> ModuleArgs:
    """ Sets up arguments for this module. Must return a ModuleArgs object.

        override_safeties is only set here to allow for the destinations of the
        arguments to not include the 'dummy' prefix for the placeholder
        arguments.
    """
    args = ModuleArgs('ClusterFinder options', 'cf')
    args.add_analysis_toggle('borders-only',
                             dest='borders_only',
                             action='store_true',
                             default=False,
                             help="Only annotate borders of existing clusters.")
    args.add_analysis_toggle('create-clusters',
                             dest='create_clusters',
                             action='store_true',
                             default=False,
                             help="Find extra clusters.")
    args.add_option('min-cds',
                    dest='min_cds_features',
                    metavar="INT",
                    type=int,
                    default=5,
                    help="Minimum size of a ClusterFinder cluster, in number of"
                         " CDS features (default: %(default)s).")
    args.add_option('mean-threshold',
                    dest='cf_threshold',
                    metavar="FLOAT",
                    type=float,
                    default=0.6,
                    help="Minimum mean probability threshold (default: %(default)s).")
    args.add_option('min-pfams',
                    dest='min_pfams',
                    metavar="INT",
                    type=int,
                    default=5,
                    help="Minimum number of biosynthetic PFam domains in a"
                         " ClusterFinder cluster (default: %(default)s).")
    return args


def is_enabled(options: ConfigType) -> bool:
    """  Uses the supplied options to determine if the module should be run
    """
    return options.cf_borders_only or options.cf_create_clusters


def check_options(options: ConfigType) -> List[str]:
    """ Checks that extra options are valid """
    if not 0. < options.cf_threshold <= 1.:
        return ["probability threshold (--cf-mean-threshold) must be between 0 and 1"]
    return []


class ClusterFinderResults(DetectionResults):
    """ Storage for predictions """
    def __init__(self, record_id: str, borders: List[ClusterBorder], create: bool = False) -> None:
        super().__init__(record_id)
        self.create_new_clusters = create
        assert isinstance(borders, list), type(borders)
        self.borders = borders

    def to_json(self) -> Dict[str, Any]:  # TODO: implement to/from json
        logging.critical("cluster_finder results always empty")
        return {}

    @staticmethod
    def from_json(json: Dict[str, Any], record: Record) -> "ClusterFinderResults":
        raise NotImplementedError("No conversion exists yet for ClusterFinderResults from JSON")

    def add_to_record(self, record: Record) -> None:
        if self.create_new_clusters:  # then get_predictions covered it already
            return
        for border in self.borders:
            record.add_cluster_border(border)

    def get_predictions(self) -> List[ClusterBorder]:
        if not self.create_new_clusters:  # then don't predict, just annotate
            return []
        return self.borders


def check_prereqs() -> List[str]:
    """Don't check for prerequisites, we don't have any"""
    return []


def run_on_record(record: Record, results: Optional[ClusterFinderResults],
                  options: ConfigType) -> ClusterFinderResults:
    """Load the data and run the cluster_predict tool"""
    if results:
        return results

    logging.info('Running ClusterFinder HMM to detect gene clusters')
    pfam_features = record.get_pfam_domains()
    if not pfam_features:
        logging.debug("No PFAM domains in record, probabilistic clusters cannot be found")
        return ClusterFinderResults(record.id, [])

    pfam_ids = []  # type: List[str]
    for pfam in pfam_features:
        pfam_ids.extend(pfam.db_xref)
    if not pfam_ids:
        logging.debug("No valid PFAM ids in record, probabilistic clusters cannot be found")
        return ClusterFinderResults(record.id, [])

    # TODO: change when PFAMs have enforced ID attributes
    pfam_features_with_ids = [feature for feature in pfam_features if feature.db_xref]

    # annotate ClusterFinder probabilities within PFAM features
    probabilities = get_pfam_probabilities(pfam_ids)
    for pfam, probabilitiy in zip(pfam_features_with_ids, probabilities):
        pfam.probability = probabilitiy

    return generate_results(record, options)


def find_rule_based_clusters(record: Record, _options: ConfigType) -> List[ClusterBorder]:
    """ Generate cluster predictions using HMM profile base rules.

        Arguments:
            record: the record to find clusters in
            options: antismash config

        Returns:
            a list of ClusterBorders, one for each matching rule
    """
    # TODO: don't go directly to hmm_detection's data
    signatures = path.get_full_path(__file__, "..", "hmm_detection", "data", "hmmdetails.txt")
    seeds = path.get_full_path(__file__, "..", "hmm_detection", "data", "bgc_seeds.hmm")
    rules = path.get_full_path(__file__, "cluster_rules.txt")
    equivalences = path.get_full_path(__file__, "..", "hmm_detection", "filterhmmdetails.txt")
    results = detect_borders_and_signatures(record, signatures, seeds, rules, equivalences,
                                            "cluster-finder")
    results.annotate_cds_features()
    logging.debug("ClusterFinder detected %d rule-based clusters", len(results.borders))
    return results.borders


def generate_results(record: Record, options: ConfigType) -> ClusterFinderResults:
    """ Find and construct cluster borders """
    rule_clusters = find_rule_based_clusters(record, options)
    prob_clusters = find_probabilistic_clusters(record, options)
    new_clusters = []
    new_clusters.extend(rule_clusters)
    for cluster in prob_clusters:
        new_cluster = ClusterBorder(cluster.location, tool="clusterfinder",
                                    probability=cluster.probability, product=PUTATIVE_PRODUCT,
                                    high_priority_product=False)
        new_clusters.append(new_cluster)
    if options.cf_create_clusters:
        for border in new_clusters:
            record.add_cluster_border(border)
    return ClusterFinderResults(record.id, new_clusters, create=options.cf_create_clusters)
