# License: GNU Affero General Public License v3 or later
# A copy of GNU AGPL v3 should have been included in this software package in LICENSE.txt.
"""
Identify conserved active site residues in PFAM_Doman / aSDomain features
"""

import logging
from typing import Any, Dict, List, Optional, Tuple

from antismash.common import module_results, path, secmet
from antismash.config import ConfigType
from antismash.config.args import ModuleArgs

from .analysis import run_all_analyses

NAME = "ActiveSiteFinder"
SHORT_DESCRIPTION = "ActiveSiteFinder identifies conserved active sites in PFAM_Domain/aSDomain features"


class ASFResults(module_results.ModuleResults):
    """ Results for active site finder """
    schema_version = 1

    def __init__(self, record_id: str, pairings: List[Tuple[secmet.features.Domain, List[str]]]) -> None:
        # pairing features will be either ModularDomain or PFAMDomain
        super().__init__(record_id)
        self.pairings = pairings

    def to_json(self) -> Dict[str, Any]:
        json = {"schema version": self.schema_version,
                "record id": self.record_id,
                "pairings": [(feature.get_name(), results) for feature, results in self.pairings]}
        return json

    @staticmethod
    def from_json(json: Dict[str, Any], record: secmet.Record) -> Optional["ASFResults"]:
        if ASFResults.schema_version != json.pop("schema version", None):
            logging.warning("Dropping ASF results, schema version has changed")
            return None
        if record.id != json.pop("record id", None):
            raise ValueError("ASF results contained mismatching record ids")

        pairings = []
        for domain_name, labels in json["pairings"]:
            domain = record.get_domain_by_name(domain_name)
            pairings.append((domain, labels))

        return ASFResults(record.id, pairings)

    def add_to_record(self, record: secmet.Record) -> None:
        assert record.id == self.record_id
        for feature, results in self.pairings:
            for result in results:
                feature.asf.add(result)


def check_options(_options: ConfigType) -> List[str]:
    """ Checks options for conflicts.
        No extra options, so they can't have conflicts.
    """
    return []


def check_prereqs(options: ConfigType) -> List[str]:
    "Checks if all required files and applications are around"
    failure_messages = []

    for binary_name in ['hmmpfam2', 'hmmscan', 'hmmpress']:
        if binary_name not in options.executables:
            failure_messages.append(f"Failed to locate file: {binary_name!r}")

    # Get all HMM profile names from XML file
    for profile in ["PKSI-KR.hmm2", "PKSI-KS_N.hmm2", "PKSI-KS_C.hmm2", "PKSI-AT.hmm2",
                    "PKSI-ACP.hmm2", "PKSI-DH.hmm2", "Thioesterase.hmm2", "PKSI-ER.hmm2",
                    "p450.hmm2"]:
        full_hmm_path = path.get_full_path(__file__, "data", profile)

        if path.locate_file(full_hmm_path) is None:
            failure_messages.append(f"Failed to locate file: {profile!r}")
            continue

    return failure_messages


def get_arguments() -> ModuleArgs:
    """ Build and return arguments. """
    args = ModuleArgs('Additional analysis', 'asf')
    args.add_analysis_toggle('--asf',
                             dest='asf',
                             action='store_true',
                             default=False,
                             help="Run active site finder analysis.")
    return args


def is_enabled(options: ConfigType) -> bool:
    """ Should the module be run with these options """
    return options.asf


def regenerate_previous_results(results: Dict[str, Any], record: secmet.Record,
                                _options: ConfigType) -> Optional[ASFResults]:
    """ Regenerate the previous results from JSON format. """
    return ASFResults.from_json(results, record)


def run_on_record(record: secmet.Record, results: Optional[ASFResults],
                  _options: ConfigType) -> ASFResults:
    """ Run the analysis, unless the previous results apply to the given record """
    if results:
        assert isinstance(results, ASFResults), type(results)
        return results
    return ASFResults(record.id, run_all_analyses(record))
