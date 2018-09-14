# License: GNU Affero General Public License v3 or later
# A copy of GNU AGPL v3 should have been included in this software package in LICENSE.txt.

""" Probability based cluster detection
"""

import logging
from typing import Any, Dict, List, Optional

from antismash.common.module_results import DetectionResults
from antismash.common.secmet import SubRegion, Record
from antismash.config import ConfigType
from antismash.config.args import ModuleArgs

from .probabilistic import find_probabilistic_clusters, get_pfam_probabilities

NAME = "clusterfinder"
SHORT_DESCRIPTION = "Statistics-based cluster detection"

PUTATIVE_PRODUCT = "cf_putative"


def get_arguments() -> ModuleArgs:
    """ Sets up arguments for this module. Must return a ModuleArgs object.
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
    if options.cf_borders_only and options.cf_create_clusters:
        return ["both cf-borders-only and cf-create-clusters specified, only one can be true"]
    return []


class ClusterFinderResults(DetectionResults):
    """ Storage for predictions """
    def __init__(self, record_id: str, areas: List[SubRegion], create: bool = False) -> None:
        super().__init__(record_id)
        self.create_new_clusters = create
        assert isinstance(areas, list), type(areas)
        self.areas = areas

    def to_json(self) -> Dict[str, Any]:  # TODO: implement to/from json
        logging.critical("cluster_finder results always empty")
        return {}

    @staticmethod
    def from_json(json: Dict[str, Any], record: Record) -> "ClusterFinderResults":
        raise NotImplementedError("No conversion exists yet for ClusterFinderResults from JSON")

    def add_to_record(self, record: Record) -> None:
        if self.create_new_clusters:  # then get_predicted_subregions covered it already
            return
        for area in self.areas:
            record.add_subregion(area)

    def get_predicted_subregions(self) -> List[SubRegion]:
        if not self.create_new_clusters:  # then don't predict, just annotate
            return []
        return self.areas


def check_prereqs() -> List[str]:
    """Don't check for prerequisites, we don't have any"""
    return []


def run_on_record(record: Record, results: Optional[ClusterFinderResults],
                  options: ConfigType) -> ClusterFinderResults:
    """Load the data and run the cluster_predict tool"""
    if results:
        return results

    logging.info('Running ClusterFinder to detect probabilistic gene clusters')
    pfam_features = record.get_pfam_domains()
    if not pfam_features:
        logging.debug("No PFAM domains in record, probabilistic clusters cannot be found")
        return ClusterFinderResults(record.id, [])

    pfam_ids = []  # type: List[str]
    for pfam in pfam_features:
        pfam_ids.append(pfam.identifier)
    if not pfam_ids:
        logging.debug("No valid PFAM ids in record, probabilistic clusters cannot be found")
        return ClusterFinderResults(record.id, [])

    # annotate ClusterFinder probabilities within PFAM features
    probabilities = get_pfam_probabilities(pfam_ids)
    for pfam, probability in zip(pfam_features, probabilities):
        pfam.probability = probability

    return generate_results(record, options)


def regenerate_previous_results(_previous: Dict[str, Any], _record: Record,
                                _options: ConfigType) -> Optional[ClusterFinderResults]:
    """ Rebuild the previous run results from a JSON object into this module's
        python results class.

        Arguments:
            previous: the previous results as a dictionary
            record: the Record that was used to generate the previous results
            options: an antismash.Config object
    """
    logging.critical("clusterfinder not regenerating results")
    return None


def generate_results(record: Record, options: ConfigType) -> ClusterFinderResults:
    """ Find and construct probabilistic cluster areas """
    predictions = find_probabilistic_clusters(record, options)
    new_areas = []
    for prediction in predictions:
        new_areas.append(SubRegion(prediction.location, tool="clusterfinder",
                                   probability=prediction.probability))
    if options.cf_create_clusters:
        for area in new_areas:
            record.add_subregion(area)
    return ClusterFinderResults(record.id, new_areas, create=options.cf_create_clusters)
