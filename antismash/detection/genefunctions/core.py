# License: GNU Affero General Public License v3 or later
# A copy of GNU AGPL v3 should have been included in this software package in LICENSE.txt.

""" Core functions and classes for detection gene functions """

import logging
from typing import Any, Dict, List, Optional

from antismash.common import fasta, subprocessing, utils
from antismash.common.hmmscan_refinement import HMMResult, refine_hmmscan_results
from antismash.common.module_results import DetectionResults
from antismash.common.secmet import CDSFeature, Record
from antismash.common.secmet.qualifiers import GeneFunction


class FunctionResults(DetectionResults):
    """ A results container for results functions for a particular tool """
    schema_version = 2

    def __init__(self, record_id: str, tool: str, best_hits: Dict[str, HMMResult],
                 function_mapping: Dict[str, GeneFunction]) -> None:
        super().__init__(record_id)
        self.tool = tool
        self.best_hits = best_hits
        self.function_mapping = function_mapping

    def add_to_record(self, record: Record) -> None:
        """ Annotate resistance genes in CDS features """
        logging.debug("annotating CDS features with %s info: %d CDSes",
                      self.tool, len(self.best_hits))
        for feature_name, result in self.best_hits.items():
            function = self.function_mapping[feature_name]
            feature = record.get_cds_by_name(feature_name)
            feature.gene_functions.add(function, self.tool,
                                       f"{result.hit_id} (Score: {result.bitscore:G}; E-value: {result.evalue})")

    @staticmethod
    def from_json(json: Dict[str, Any], record: Record) -> Optional["FunctionResults"]:
        if json.get("schema_version") != FunctionResults.schema_version:
            logging.debug("Schema version mismatch, discarding FunctionResults for tool: %s",
                          json.get("tool", "unknown"))
            return None
        if record.id != json.get("record_id"):
            logging.debug("Record ID mismatch, discarding FunctionResults for tool: %s",
                          json.get("tool", "unknown"))
            return None
        hits = {}
        for hit, parts in json["best_hits"].items():
            hits[hit] = HMMResult.from_json(parts)

        mapping = {}
        for cds_name, simple_function in json["mapping"].items():
            mapping[cds_name] = GeneFunction.from_string(simple_function)

        results = FunctionResults(json["record_id"], json["tool"], hits, mapping)
        return results

    def to_json(self) -> Dict[str, Any]:
        return {"schema_version": self.schema_version,
                "record_id": self.record_id,
                "tool": self.tool,
                "best_hits": {key: val.to_json() for key, val in self.best_hits.items()},
                "mapping": {name: str(function) for name, function in self.function_mapping.items()},
                }


def scan_for_functions(cds_features: List[CDSFeature], database: str,
                       hmmscan_opts: Optional[List[str]] = None) -> Dict[str, HMMResult]:
    """ Finds possible classifications for the provided genes.

        Arguments:
            cds_features: a list of CDSFeatures to classify
            database: the path to the database to check
            hmmscan_opts: a list of extra options to provide to hmmscan

        Returns:
            a dictionary mapping CDS name to a list of HMMResult instances of
                classifications
    """
    search_fasta = fasta.get_fasta_from_features(cds_features)
    results = subprocessing.run_hmmscan(database, search_fasta, hmmscan_opts)
    hmm_lengths = utils.get_hmm_lengths(database)
    hmm_results = refine_hmmscan_results(results, hmm_lengths)

    best_hits: Dict[str, HMMResult] = {}

    for cds in cds_features:
        cds_name = cds.get_name()
        hits = hmm_results.get(cds_name)
        if not hits:
            continue
        best_hits[cds_name] = hits[0]

    return best_hits
