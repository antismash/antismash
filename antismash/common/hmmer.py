# License: GNU Affero General Public License v3 or later
# A copy of GNU AGPL v3 should have been included in this software package in LICENSE.txt.

""" Finds all PFAM domains found in features throughout a record. """

import logging
import os
from typing import Any, Dict, List, Optional, Tuple

from antismash.common import fasta, module_results, pfamdb, subprocessing
from antismash.common.secmet import Record, CDSFeature
from antismash.common.secmet.feature import FeatureLocation, PFAMDomain


class HmmerResults(module_results.ModuleResults):
    """ Results for full-hmmer """
    schema_version = 1

    def __init__(self, record_id: str, evalue: float, score: float,
                 database: str, tool: str, hits: List[Dict[str, Any]]) -> None:
        super().__init__(record_id)
        self.hits = list(hits)
        self.evalue = float(evalue)
        self.score = float(score)
        self.database = str(database)
        self.tool = str(tool)

    def to_json(self) -> Dict[str, Any]:
        json = {"hits": self.hits, "record id": self.record_id,
                "schema": self.schema_version, "max evalue": self.evalue,
                "min score": self.score, "database": self.database,
                "tool": self.tool}
        return json

    @staticmethod
    def from_json(json: Dict[str, Any], record: Record, max_evalue: float,
                  min_score: float, database: str) -> Optional["HmmerResults"]:
        """ Regenerate the results from JSON.
            If max_evalue or min_score aren't equal or narrower than those the
            results were generated with, the results will be discarded.
        """

        if record.id != json.get("record id"):
            logging.warning("Hmmer results are for different record, discarding previous results")
            return None

        if json.get("schema") != HmmerResults.schema_version:
            logging.warning("Hmmer results are for different result schema, discarding previous results")
            return None

        if json.get("database") != database:
            logging.warning("Hmmer database has changed, discarding previous results")
            return None

        evalue = json.get("max evalue")
        score = json.get("min score")
        if evalue is None or score is None:
            raise ValueError("Invalid Hmmer result values")
        assert isinstance(score, float) and isinstance(evalue, float)

        # if the current options have expanded the detection range, rerunning is required
        if evalue < max_evalue:
            logging.debug("Discarding Hmmer results, using new evalue threshold")
            return None
        if score > min_score:
            logging.debug("Discarding Hmmer results, using new score threshold")
            return None

        hits = json.get("hits")
        if not isinstance(hits, list):
            raise TypeError("FullHmmer results contain unexpected types")
        # if the thresholds changed, trim out any extra hits here
        hits = [hit for hit in hits if hit["score"] >= min_score and hit["evalue"] <= max_evalue]

        results = HmmerResults(record.id, max_evalue, min_score, json["database"], json["tool"], hits)
        return results

    def add_to_record(self, record: Record) -> None:
        db_version = pfamdb.get_db_version_from_path(self.database)
        for i, hit in enumerate(self.hits):
            pfam_feature = PFAMDomain(FeatureLocation(hit["start"], hit["end"], hit["strand"]),
                                      description=hit["description"])
            for key in ["label", "locus_tag", "domain", "evalue",
                        "score", "translation", "db_xref"]:
                setattr(pfam_feature, key, hit[key])
            pfam_feature.tool = self.tool
            pfam_feature.database = db_version
            pfam_feature.detection = "hmmscan"
            pfam_feature.domain_id = "{}_{}_{:04d}".format(self.tool, pfam_feature.locus_tag, i + 1)
            record.add_pfam_domain(pfam_feature)


def calculate_start_and_end(feature, result) -> Tuple[int, int]:
    "Calculate start and end of a result"
    if feature.strand == 1:
        start = feature.location.start + (3 * result.query_start)
        end = feature.location.start + (3 * result.query_end)
    else:
        end = feature.location.end - (3 * result.query_start)
        start = feature.location.end - (3 * result.query_end)

    return start, end


def build_hits(record, hmmscan_results, min_score: float, max_evalue: float, database: str) -> List[Dict[str, Any]]:
    "Builds PFAMDomains from the given hmmscan results"
    logging.debug("Generating feature objects for PFAM hits")

    hits = []
    feature_by_id = record.get_cds_name_mapping()

    for result in hmmscan_results:
        for hsp in result.hsps:
            if hsp.bitscore <= min_score or hsp.evalue >= max_evalue:
                continue

            if hsp.query_id not in hsp.query_id:
                continue

            feature = feature_by_id[hsp.query_id]

            start, end = calculate_start_and_end(feature, hsp)

            dummy_feature = PFAMDomain(FeatureLocation(start, end, feature.location.strand),
                                       description="")

            hit = {"start": start, "end": end, "strand": feature.location.strand,
                   "label": result.id, "locus_tag": feature.locus_tag,
                   "domain": hsp.hit_id, "evalue": hsp.evalue, "score": hsp.bitscore,
                   "translation": str(dummy_feature.extract(record.seq).translate(table=feature.transl_table)),
                   "db_xref": [pfamdb.get_pfam_id_from_name(hsp.hit_id, database)],
                   "description": hsp.hit_description}
            hits.append(hit)
    return hits


def run_hmmer(record: Record, features: List[CDSFeature], max_evalue: float,
              min_score: float, database: str, tool: str) -> HmmerResults:
    """ Build hmmer results for the given features"""
    if not os.path.exists(database):
        raise ValueError("Given database does not exist: %s" % database)
    query_sequence = fasta.get_fasta_from_features(features)
    hmmscan_results = subprocessing.run_hmmscan(database, query_sequence)
    hits = build_hits(record, hmmscan_results, min_score, max_evalue, database)
    return HmmerResults(record.id, max_evalue, min_score, database, tool, hits)
