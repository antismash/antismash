# License: GNU Affero General Public License v3 or later
# A copy of GNU AGPL v3 should have been included in this software package in LICENSE.txt.

""" Common functionality for finding and marking PFAMDomains within a record. """

import logging
import os
from typing import Any, Dict, Iterable, List, Optional

from antismash.common import fasta, module_results, path, pfamdb, subprocessing
from antismash.common.secmet import Record, CDSFeature
from antismash.common.secmet.features import PFAMDomain
from antismash.common.secmet.locations import location_from_string


class HmmerResults(module_results.ModuleResults):
    """ Results for hmmer-based detection """
    schema_version = 2

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
    def from_json(json: Dict[str, Any], record: Record, max_evalue: float,  # type: ignore  # pylint: disable=arguments-differ
                  min_score: float) -> Optional["HmmerResults"]:
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
            pfam_feature = PFAMDomain(location_from_string(hit["location"]),
                                      description=hit["description"], protein_start=hit["protein_start"],
                                      protein_end=hit["protein_end"], identifier=hit["identifier"],
                                      tool=self.tool)
            for key in ["label", "locus_tag", "domain", "evalue",
                        "score", "translation"]:
                setattr(pfam_feature, key, hit[key])
            pfam_feature.database = db_version
            pfam_feature.detection = "hmmscan"
            pfam_feature.domain_id = "{}_{}_{:04d}".format(self.tool, pfam_feature.locus_tag, i + 1)
            record.add_pfam_domain(pfam_feature)


def build_hits(record: Record, hmmscan_results: List, min_score: float,
               max_evalue: float, database: str) -> List[Dict[str, Any]]:
    """ Builds PFAMDomains from the given hmmscan results

        Arguments:
            record: the Record being scanned
            hmmscan_results: the results of Bio.SearchIO.parse
            min_score: a minimum allowable bitscore for hits (exclusive)
            max_evalue: a maximum allowable evalue for hits (exclusive)
            database: the name of the database used to find the hits

        Returns:
            a list of JSON representations of hmmer hits
    """
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
            location = feature.get_sub_location_from_protein_coordinates(hsp.query_start, hsp.query_end)

            hit = {"location": str(location),
                   "label": result.id, "locus_tag": feature.get_name(),
                   "domain": hsp.hit_id, "evalue": hsp.evalue, "score": hsp.bitscore,
                   "translation": feature.translation[hsp.query_start:hsp.query_end + 1],
                   "identifier": pfamdb.get_pfam_id_from_name(hsp.hit_id, database),
                   "description": hsp.hit_description, "protein_start": hsp.query_start, "protein_end": hsp.query_end}
            hits.append(hit)
    return hits


def run_hmmer(record: Record, features: Iterable[CDSFeature], max_evalue: float,
              min_score: float, database: str, tool: str) -> HmmerResults:
    """ Build hmmer results for the given features

        Arguments:
            record: the Record instance to run hmmer over
            features: the list of CDSFeatures to run over specifically
            max_evalue: a maximum evalue allowed for hits (exclusive)
            min_evalue: a minimum evalue allowed for hits (exclusive)
            database: the database to search for hits within
            tool: the name of the specific tool calling into this module
    """
    if not os.path.exists(database):
        raise ValueError("Given database does not exist: %s" % database)
    query_sequence = fasta.get_fasta_from_features(features)
    hmmscan_results = subprocessing.run_hmmscan(database, query_sequence, opts=["--cut_tc"])
    hits = build_hits(record, hmmscan_results, min_score, max_evalue, database)
    return HmmerResults(record.id, max_evalue, min_score, database, tool, hits)


def ensure_database_pressed(filepath: str, return_not_raise: bool = False) -> List[str]:
    """ Ensures that the given HMMer database exists and that the hmmpress
        generated files aren't out of date.

        Arguments:
            filepath: the path to the HMMer database
            return_not_raise: whether to catch errors and return their messages as strings

        Returns:
            any encountered error messages, will never be populated without return_not_raise == True
    """
    try:
        modified_time = os.path.getmtime(filepath)
    except FileNotFoundError as err:
        if not return_not_raise:
            raise
        return [str(err)]
    components = ["{}{}".format(filepath, ext) for ext in ['.h3f', '.h3i', '.h3m', '.h3p']]
    outdated = False
    for component in components:
        if not path.locate_file(component) or os.path.getmtime(component) < modified_time:
            logging.info("%s does not exist or is out of date, hmmpressing %s",
                         component, filepath)
            outdated = True
            break

    if outdated:
        result = subprocessing.run_hmmpress(filepath)
        if not result.successful():
            msg = "Failed to hmmpress {!r}: {}".format(filepath, result.stderr)
            if not return_not_raise:
                raise RuntimeError(msg)
            return [msg]
    return []
