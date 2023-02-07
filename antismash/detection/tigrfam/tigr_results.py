# License: GNU Affero General Public License v3 or later
# A copy of GNU AGPL v3 should have been included in this software package in LICENSE.txt.

""" Results class for the TIGRFam detection module. """

from typing import Any, Dict, List, Optional

from antismash.common.hmmer import HmmerResults, HmmerHit
from antismash.common.secmet import Record
from antismash.common.secmet.features import FeatureLocation
from antismash.common.secmet.locations import location_from_string

from .tigr_domain import TIGRDomain


class TIGRFamResults(HmmerResults):
    """ Results class for the TIGRFam detecion module, similar to HmmerResults, but not bound to PFAM"""
    def __init__(self, record_id: str, evalue: float, score: float, hits: List[HmmerHit]) -> None:
        super().__init__(record_id, evalue, score, "TIGRFam.hmm", "TIGRFam", hits)

    def add_to_record(self, record: Record) -> None:
        """ Adds the hits as TIGRDomains to the given record """
        if record.id != self.record_id:
            raise ValueError("Record to store in and record analysed don't match")

        for i, hit in enumerate(self.hits):
            protein_location = FeatureLocation(hit.protein_start, hit.protein_end)
            tigr_feature = TIGRDomain(location_from_string(hit.location),
                                      description=hit.description, protein_location=protein_location,
                                      identifier=hit.identifier, locus_tag=hit.locus_tag)
            for key in ["label", "locus_tag", "domain", "evalue",
                        "score", "translation"]:
                setattr(tigr_feature, key, getattr(hit, key))
            tigr_feature.detection = "hmmscan"
            tigr_feature.domain_id = f"{self.tool}_{tigr_feature.locus_tag}_{i + 1:04d}"
            record.add_feature(tigr_feature)

    def refilter(self, max_evalue: float, min_score: float) -> "TIGRFamResults":
        return TIGRFamResults.from_hmmer_results(super().refilter(max_evalue, min_score))

    @classmethod
    def from_json(cls, json: Dict[str, Any], record: Record) -> Optional["TIGRFamResults"]:
        results = super().from_json(json, record)
        if results is None:
            return None
        return cls.from_hmmer_results(results)

    @classmethod
    def from_hmmer_results(cls, hmmer_results: HmmerResults) -> "TIGRFamResults":
        """ Convert HmmerResults into TIGRFamResults. """
        return cls(hmmer_results.record_id, hmmer_results.evalue, hmmer_results.score, hmmer_results.hits)
