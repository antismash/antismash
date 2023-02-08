# License: GNU Affero General Public License v3 or later
# A copy of GNU AGPL v3 should have been included in this software package in LICENSE.txt.

"""Integration of precision mode of RREFinder
   for the detection of RiPP Recognition Elements
   (RREs) in RiPP BGCs
"""

import logging

from collections import defaultdict
from typing import Any, Dict, List, Optional, Tuple

from antismash.config import get_config
from antismash.common.hmmer import run_hmmer, HmmerResults, HmmerHit
from antismash.common.module_results import ModuleResults
from antismash.common.secmet import Record
from antismash.common.secmet.features import CDSFeature, Feature, FeatureLocation
from antismash.common.secmet.locations import location_from_string

from .rre_domain import RREDomain


class RREFinderResults(ModuleResults):
    """ Results class for the RREFinder analysis"""
    schema_version = 1

    def __init__(self, record_id: str, bitscore_cutoff: float, min_length: int,
                 hits_by_protocluster: Dict[int, List[str]], hits_by_cds: Dict[str, List[HmmerHit]]
                 ) -> None:
        super().__init__(record_id)
        # The cutoffs used for hmmscan
        self.bitscore_cutoff = bitscore_cutoff
        self.min_length = min_length
        # All the locus tags that are an RRE hit per protocluster
        self.hits_by_protocluster = hits_by_protocluster
        # All the hit info per locus tag
        self.hits_by_cds = hits_by_cds
        self.features: List[Feature] = []  # features created for RREs
        self.tool = 'rrefinder'
        self.database = 'RREFam.hmm'
        self.detection = 'hmmscan'
        self.convert_hits_to_features()

    def convert_hits_to_features(self) -> None:
        """Convert all the hits found to features"""
        for locus_tag, hits in self.hits_by_cds.items():
            domain_counts: Dict[str, int] = defaultdict(int)
            for hit in hits:
                location = location_from_string(hit.location)
                protein_location = FeatureLocation(hit.protein_start, hit.protein_end)
                rre_feature = RREDomain(location, hit.description, protein_location,
                                        identifier=hit.identifier, locus_tag=locus_tag, domain=hit.domain)

                # Set additional properties
                rre_feature.score = hit.score
                rre_feature.evalue = hit.evalue
                rre_feature.label = hit.label
                rre_feature.translation = hit.translation

                rre_feature.database = self.database
                rre_feature.detection = self.detection

                domain_counts[hit.domain] += 1  # 1-indexed, so increment before use
                rre_feature.domain_id = f"{self.tool}_{locus_tag}_{hit.domain}.{domain_counts[hit.domain]}"

                self.features.append(rre_feature)

    def add_to_record(self, record: Record) -> None:
        """ Adds the analysis results to the record """
        if record.id != self.record_id:
            raise ValueError("Record to store in and record analysed don't match")

        for feature in self.features:
            record.add_feature(feature)

    def to_json(self) -> Dict[str, Any]:
        """ Constructs a JSON representation of this instance """

        by_cds = {}
        for name, hits in self.hits_by_cds.items():
            by_cds[name] = [hit.to_json() for hit in hits]

        return {
            "schema_version": self.schema_version,
            "bitscore_cutoff": self.bitscore_cutoff,
            "hits_by_protocluster": self.hits_by_protocluster,
            "hits_by_cds": by_cds,
            "min_length": self.min_length,
            "record_id": self.record_id,
        }

    @staticmethod
    def from_json(json: Dict[str, Any], record: Record) -> Optional["RREFinderResults"]:
        """ Regenerates the RREFinderResults from json.
            If less strict settings are given (e.g. a lower cutoff
            or a shorter minimum length), the results are discarded.
        """

        # check that the previous data version is the same as current, if not, discard the results
        if json["schema_version"] != RREFinderResults.schema_version:
            return None

        if record.id != json.get("record_id"):
            logging.warning("RREFinder results are for different record, discarding previous results")
            return None

        prev_min_length = json["min_length"]
        prev_bitscore_cutoff = json["bitscore_cutoff"]

        assert isinstance(prev_min_length, int) and isinstance(prev_bitscore_cutoff, float)

        # Grab the current cutoffs
        options = get_config()
        bitscore_cutoff = options.rre_cutoff
        min_length = options.rre_min_length

        # Check if the cutoff options set in this run are the same or stricter than in the previous run
        if bitscore_cutoff < prev_bitscore_cutoff:
            logging.debug("RREFinderResults bitscore cutoff has changed, discarding previous results")
            return None
        if min_length < prev_min_length:
            logging.debug("RREFinderResults minimum length has changed, discarding previous results")
            return None

        hits_by_cds = {}
        for name, hits in json["hits_by_cds"].items():
            hits_by_cds[name] = [HmmerHit.from_json(hit) for hit in hits]

        # Refilter the hits (in case the cutoff is now more stringent)
        by_proto = {int(key): val for key, val in json['hits_by_protocluster'].items()}
        filtered_hits_by_cds, filtered_hits_by_protocluster = filter_hits(hits_by_cds, by_proto,
                                                                          min_length, bitscore_cutoff)
        return RREFinderResults(record.id, bitscore_cutoff, min_length,
                                filtered_hits_by_protocluster, filtered_hits_by_cds)

    def refilter(self, min_length: int, min_score: float) -> "RREFinderResults":
        """Trims the results to stricter thresholds for score and length"""
        if min_length < self.min_length:
            raise ValueError(f"cannot refilter to a more lenient length: {self.min_length} -> {min_length}")
        if min_score < self.bitscore_cutoff:
            raise ValueError(f"cannot refilter to a more lenient score: {self.bitscore_cutoff} -> {min_score}")
        by_cds, by_proto = filter_hits(self.hits_by_cds, self.hits_by_protocluster, min_length, min_score)
        self.hits_by_cds = by_cds
        self.hits_by_protocluster = by_proto
        return self


def gather_rre_candidates(record: Record) -> Tuple[Dict[int, List[str]], Dict[str, CDSFeature]]:
    '''Gather all RRE candidates that need to be analyzed with hmmscan
       and all unique candidates (by CDS name) to prevent double analysis
       and features in the case of overlapping RiPP protoclusters.
    '''
    rre_candidates_by_protocluster: Dict[int, List[str]] = defaultdict(list)
    cds_info: Dict[str, CDSFeature] = {}

    for region in record.get_regions():
        for protocluster in region.get_unique_protoclusters():
            if protocluster.product_category == "RiPP":
                protocluster_number = protocluster.get_protocluster_number()
                for cds in protocluster.cds_children:
                    cds_name = cds.get_name()
                    rre_candidates_by_protocluster[protocluster_number].append(cds_name)
                    cds_info[cds_name] = cds
    return dict(rre_candidates_by_protocluster), cds_info


def extract_rre_hits(hmm_result: HmmerResults) -> Dict[str, List[HmmerHit]]:
    '''Extract the hits per locus_tag from a HmmerResults object'''
    hits_by_cds: Dict[str, List[HmmerHit]] = defaultdict(list)
    for hit in hmm_result.hits:
        hits_by_cds[hit.locus_tag].append(hit)
    return dict(hits_by_cds)


def filter_hits(hits_by_cds: Dict[str, List[HmmerHit]], candidates_by_protocluster: Dict[int, List[str]],
                min_length: int, bitscore_cutoff: float
                ) -> Tuple[Dict[str, List[HmmerHit]], Dict[int, List[str]]]:
    '''Filter the hits based on the bitscore and length criteria'''
    filtered_hits_by_cds: Dict[str, List[HmmerHit]] = {}
    for cds_name, hits in hits_by_cds.items():
        trimmed_hits = [hit for hit in hits if check_hmm_hit(hit, min_length, bitscore_cutoff)]
        if trimmed_hits:
            filtered_hits_by_cds[cds_name] = trimmed_hits

    filtered_tags_by_protocluster: Dict[int, List[str]] = {}
    for protocluster, locus_tags in candidates_by_protocluster.items():
        trimmed = [locus for locus in locus_tags if locus in filtered_hits_by_cds]
        if trimmed:
            filtered_tags_by_protocluster[protocluster] = trimmed
    return filtered_hits_by_cds, filtered_tags_by_protocluster


def check_hmm_hit(hit: HmmerHit, min_length: int, bitscore_cutoff: float) -> bool:
    """Check if a HMM hit passes the set cutoffs
    """
    return hit.score >= bitscore_cutoff and len(hit) >= min_length


def run_rrefinder(record: Record, bitscore_cutoff: float, min_length: int, database: str) -> RREFinderResults:
    """Run RREFinder on a given record
    """
    # Gather all RRE candidates
    candidates_by_protocluster, cds_info = gather_rre_candidates(record)
    # Run hmmscan per protocluster and gather the hits
    if not cds_info:
        filtered_hits_by_protocluster: Dict[int, List[str]] = {}
        filtered_hits_by_cds: Dict[str, List[HmmerHit]] = {}
    else:
        hmm_results = run_hmmer(record, cds_info.values(), max_evalue=1, min_score=bitscore_cutoff,
                                database=database, tool='rrefinder', use_cut_tc=False,
                                filter_overlapping=False)
        # Extract the RRE hits
        hits_by_cds = extract_rre_hits(hmm_results)
        # Filter the hits
        filtered_hits_by_cds, filtered_hits_by_protocluster = filter_hits(hits_by_cds, candidates_by_protocluster,
                                                                          min_length, bitscore_cutoff)
    return RREFinderResults(record.id, bitscore_cutoff, min_length,
                            filtered_hits_by_protocluster, filtered_hits_by_cds)
