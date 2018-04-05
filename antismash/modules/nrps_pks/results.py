# License: GNU Affero General Public License v3 or later
# A copy of GNU AGPL v3 should have been included in this software package in LICENSE.txt.

""" Contains the results classes for the nrps_pks module """

import logging
from typing import Any, Dict

from antismash.common.module_results import ModuleResults

from .at_analysis.at_analysis import ATSignatureResults
from .parsers import LONG_TO_SHORT

DOMAIN_TYPE_MAPPING = {'Condensation_DCL': 'Condensation',
                       'Condensation_LCL': 'Condensation',
                       'Condensation_Dual': 'Condensation',
                       'Condensation_Starter': 'Condensation',
                       'CXglyc': 'Condensation',
                       'Cglyc': 'Condensation',
                       'cMT': 'MT',
                       'oMT': 'MT',
                       'nMT': 'MT',
                       'Polyketide_cyc': 'Polyketide_cyc',
                       'Polyketide_cyc2': 'Polyketide_cyc'}


_UNKNOWN = "(unknown)"


class PKSResults:
    """ Results for the PKS section of the nrps_pks module """
    __slots__ = ["method_results"]

    class EnforcedDict(dict):
        """ Enforces the value of the key 'signature' to be a consistent type """
        def __setitem__(self, key, val):
            if key == "signature":
                assert isinstance(val, ATSignatureResults)
            return super().__setitem__(key, val)

    def __init__(self):
        self.method_results = PKSResults.EnforcedDict()

    def to_json(self) -> Dict[str, Any]:
        """ Store results as a JSON object """
        results = {}
        for key, value in self.method_results.items():
            if key == "signature":
                value = value.to_json()
            results[key] = value
        return results

    @staticmethod
    def from_json(json: Dict[str, Any]) -> "PKSResults":
        """ Reconstruct a PKSResults instance from JSON object """
        assert isinstance(json, dict)
        results = PKSResults()
        results.method_results = json
        results.method_results["signature"] = ATSignatureResults.from_json(json.get("signature", {}))
        return results


class NRPS_PKS_Results(ModuleResults):
    """ The combined results of the nrps_pks module """
    _schema_version = 1
    __slots__ = ["_pks", "nrps", "consensus", "consensus_transat", "cluster_predictions"]

    def __init__(self, record_id):
        super().__init__(record_id)
        self._pks = PKSResults()
        self.nrps = {}  # not currently used, but will be
        self.consensus = {}
        self.cluster_predictions = {}  # type: Dict[int, str]
        self.consensus_transat = {}

    @property
    def pks(self) -> PKSResults:
        """ The result of the PKS analyses """
        return self._pks

    @pks.setter
    def pks(self, value: PKSResults):
        assert isinstance(value, PKSResults)
        self._pks = value

    def to_json(self) -> Dict[str, Any]:
        results = {"schema_version": self._schema_version}
        results["record_id"] = self.record_id
        results["pks"] = self.pks.to_json()
        results["nrps"] = self.nrps
        results["consensus"] = self.consensus
        results["cluster_predictions"] = self.cluster_predictions
        return results

    @staticmethod
    def from_json(json: Dict[str, Any], _record) -> "NRPS_PKS_Results":
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

    def add_to_record(self, record) -> None:
        """ Save substrate specificity predictions in NRPS/PKS domain sec_met info of record
        """

        for cds_feature in record.get_nrps_pks_cds_features():
            nrps_qualifier = cds_feature.nrps_pks
            for domain in nrps_qualifier.domains:
                feature = record.get_domain_by_name(domain.feature_name)
                domain.predictions.clear()
                if domain.name == "AMP-binding":
                    # no NRPS predictors right now, so all A domains are 'nrp'
                    domain.predictions["consensus"] = "nrp"

                elif domain.name == "PKS_AT":
                    # For t1pks, t2pks and t3pks
                    if 'transatpks' not in cds_feature.cluster.products:
                        consensus = self.consensus[domain.feature_name]
                    else:  # for transatpks
                        consensus = self.consensus_transat[domain.feature_name]
                    pks_sig = self.pks.method_results["signature"][domain.feature_name]
                    if pks_sig:
                        domain.predictions["PKS signature"] = pks_sig[0].name.rsplit("_", 1)[1]
                    else:
                        domain.predictions["PKS signature"] = _UNKNOWN
                    minowa = self.pks.method_results["minowa_at"][domain.feature_name][0][0]
                    domain.predictions["Minowa"] = LONG_TO_SHORT.get(minowa, minowa)
                    domain.predictions["consensus"] = consensus

                elif domain.name == "CAL_domain":
                    minowa = self.pks.method_results["minowa_cal"][domain.feature_name][0][0]
                    domain.predictions["Minowa"] = LONG_TO_SHORT.get(minowa, minowa)

                elif domain.name == "PKS_KR":
                    domain.predictions["KR activity"] = \
                            "active" if self.pks.method_results["kr_activity"][domain.feature_name] else "inactive"
                    domain.predictions["KR stereochemistry"] = \
                            self.pks.method_results["kr_stereochem"].get(domain.feature_name, _UNKNOWN)
                for method, pred in domain.predictions.items():
                    feature.specificity.append("%s: %s" % (method, pred))

                mapping = DOMAIN_TYPE_MAPPING.get(domain.name)
                if mapping:
                    feature.domain_subtype = domain.name
                    feature.domain = mapping
