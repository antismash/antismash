# License: GNU Affero General Public License v3 or later
# A copy of GNU AGPL v3 should have been included in this software package in LICENSE.txt.

import logging

from antismash.common.module_results import ModuleResults

from .at_analysis.at_analysis import ATSignatureResults

class NRPS_PKS_Results(ModuleResults):
    _schema_version = 1
    __slots__ = ["_pks", "nrps", "consensus", "consensus_transat", "cluster_predictions"]

    def __init__(self, record_id):
        super().__init__(record_id)
        self._pks = PKSResults()
        self.nrps = {}  # not currently used, but will be
        self.consensus = {}
        self.cluster_predictions = {}  # type: Dict[int, str]

    @property
    def pks(self):
        return self._pks

    @pks.setter
    def pks(self, value):
        assert isinstance(value, PKSResults)
        self._pks = value

    def to_json(self):
        results = {"schema_version": self._schema_version}
        results["record_id"] = self.record_id
        results["pks"] = self.pks.to_json()
        results["nrps"] = self.nrps
        results["consensus"] = self.consensus
        results["cluster_predictions"] = self.cluster_predictions
        return results

    @staticmethod
    def from_json(json):
        assert "record_id" in json
        if json.get("schema_version") != NRPS_PKS_Results._schema_version:
            logging.warning("Mismatching schema version, dropping results")
            return None
        results = NRPS_PKS_Results(json["record_id"])
        pks = json.get("pks")
        if pks:
            results.pks = PKSResults.from_json(pks)
        results.nrps = json.get("nrps", {})
        results.consensus = json["consensus"]
        results.cluster_predictions = json["cluster_predictions"]
        return results

    def add_to_record(self, record):
        logging.critical("NRPS/PKS results skipping add_to_record")
        return {}


class PKSResults:
    __slots__ = ["method_results"]

    class EnforcedDict(dict):
        def __setitem__(self, key, val):
            if key == "signature":
                assert isinstance(val, ATSignatureResults)
            return super().__setitem__(key, val)

    def __init__(self):
        self.method_results = PKSResults.EnforcedDict()

    def to_json(self):
        results = {}
        for key, value in self.method_results.items():
            if key == "signature":
                value = value.to_json()
            results[key] = value
        return results

    @staticmethod
    def from_json(json):
        assert isinstance(json, dict)
        results = PKSResults()
        results.method_results = json
        results.method_results["signature"] = ATSignatureResults.from_json(json.get("signature", {}))
        return results
