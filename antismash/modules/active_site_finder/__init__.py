# License: GNU Affero General Public License v3 or later
# A copy of GNU AGPL v3 should have been included in this software package in LICENSE.txt.
"""
Identify conserved active site residues in PFAM_Doman / aSDomain features
"""

import datetime
import glob
import logging
import os
from typing import Any, Dict, List, Tuple

from antismash.common import path, subprocessing, module_results, secmet
from antismash.config.args import ModuleArgs

from .analysis import run_all_analyses

NAME = "ActiveSiteFinder"
SHORT_DESCRIPTION = "ActiveSiteFinder identifies conserved active sites in PFAM_Domain/aSDomain features"


class ASFResults(module_results.ModuleResults):
    schema_version = 1
    def __init__(self, record_id: str, pairings: List[Tuple[secmet.feature.AntismashFeature, List[str]]]) -> None:
        # pairing features will be either AntismashDomain or PFAMDomain
        super().__init__(record_id)
        self.pairings = pairings

    def to_json(self) -> Dict[str, List[str]]:
        json = {"schema version": self.schema_version,
                "record id": self.record_id}
        for feature, results in self.pairings:
            json[feature.domain_id] = (feature.type, results)
        return json

    @staticmethod
    def from_json(json: Dict[str, Any], record: secmet.Record) -> "ASFResults":
        if ASFResults.schema_version != json.pop("schema version", None):
            logging.warning("Dropping ASF results, schema version has changed")
            return None
        if record.id != json.pop("record id", None):
            raise ValueError("ASF results contained mismatching record ids")
        # TODO: fetch the domain from the record using domain id
        return ASFResults(record.id, list(json.values()))

    def add_to_record(self, record):
        assert record.id == self.record_id
        for feature, results in self.pairings:
            for result in results:
                feature.asf.add(result)


def check_options(options) -> List[str]:
    return []  # TODO: maybe bail if no full_hmmer/cluster_hmmer?


def check_prereqs() -> List[str]:
    "Checks if all required files and applications are around"
    _binary_extensions = ['.h3f', '.h3i', '.h3m', '.h3p']

    failure_messages = []

    for binary_name in ['hmmpfam2', 'hmmscan', 'hmmpress']:
        if not path.locate_executable(binary_name):
            failure_messages.append("Failed to locate file: %r" % binary_name)

    # Get all HMM profile names from XML file
    for profile in ["PKSI-KR.hmm2", "PKSI-KS_N.hmm2", "PKSI-KS_C.hmm2", "PKSI-AT.hmm2",
                    "PKSI-ACP.hmm2", "PKSI-DH.hmm2", "Thioesterase.hmm2", "PKSI-ER.hmm2",
                    "aa-activating.aroundLys.hmm2", "p450.hmm2"]:
        full_hmm_path = path.get_full_path(__file__, "data", profile)

        if path.locate_file(full_hmm_path) is None:
            failure_messages.append("Failed to locate file: %s" % profile)
            continue

        if profile.endswith(".hmm2"):
            continue

        for ext in _binary_extensions:
            binary = "{hmm}{ext}".format(hmm=full_hmm_path, ext=ext)
            if not path.locate_file(binary):
                result = subprocessing.run_hmmpress(full_hmm_path)
                if not result.successful():
                    failure_messages.append("Failed to hmmpress {!r}: {!r}".format(profile, result.stderr))

                # hmmpress generates _all_ binary files in one go, so stop the loop
                break

            binary_mtime = os.path.getmtime(binary)
            hmm_mtime = os.path.getmtime(full_hmm_path)
            if hmm_mtime < binary_mtime:
                # generated file younger than hmm profile, do nothing
                continue
            try:
                for filename in glob.glob("{}.h3?".format(full_hmm_path)):
                    logging.debug("removing outdated file %r", filename)
                    os.remove(filename)
            except OSError as err:
                failure_messages.append("Failed to remove outdated binary file for %s: %s" %
                                        (profile, err))
                break
            result = subprocessing.run_hmmpress(full_hmm_path)
            if not result.successful():
                failure_messages.append("Failed to hmmpress %r: %r" % (profile, result.stderr))
                failure_messages.append("HMM binary files outdated. %s (changed: %s) vs %s (changed: %s)" %
                                        (profile, datetime.datetime.fromtimestamp(hmm_mtime),
                                         binary, datetime.datetime.fromtimestamp(binary_mtime)))
            # hmmpress generates _all_ binary files in one go, so stop the loop
            break

    return failure_messages


def get_arguments():
    args = ModuleArgs('Additional analysis', 'asf')
    args.add_analysis_toggle('--asf',
                             dest='asf',
                             action='store_true',
                             default=False,
                             help="Run active site finder analysis.")
    return args


def is_enabled(options):
    return options.asf


def regenerate_previous_results(results, record, _options) -> ASFResults:
    return ASFResults.from_json(results, record)


def run_on_record(record, results, _options) -> ASFResults:
    if results:
        assert isinstance(results, ASFResults), type(results)
        return results
    return ASFResults(record.id, run_all_analyses(record))
