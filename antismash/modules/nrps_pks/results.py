# License: GNU Affero General Public License v3 or later
# A copy of GNU AGPL v3 should have been included in this software package in LICENSE.txt.

""" Contains the results classes for the nrps_pks module """

from collections import defaultdict
import logging
from typing import Any, Dict, Tuple, Union  # pylint: disable=unused-import

from antismash.common.module_results import ModuleResults
from antismash.common.secmet import Record, AntismashDomain

from .at_analysis.at_analysis import ATSignatureResults
from .parsers import LONG_TO_SHORT, generate_nrps_consensus
from .data_structures import Prediction
from .nrps_predictor import PredictorSVMResult

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
        def __setitem__(self, key: str, val: Union[Dict[str, str], Dict[str, bool]]) -> None:
            if key == "signature":
                assert isinstance(val, ATSignatureResults)
            return super().__setitem__(key, val)

    def __init__(self) -> None:
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
        results.method_results.update(json)
        results.method_results["signature"] = ATSignatureResults.from_json(json.get("signature", {}))
        return results


class NRPS_PKS_Results(ModuleResults):
    """ The combined results of the nrps_pks module """
    _schema_version = 1
    __slots__ = ["_pks", "_nrps", "consensus", "consensus_transat", "cluster_predictions", "domain_predictions"]

    def __init__(self, record_id: str) -> None:
        super().__init__(record_id)
        self._pks = PKSResults()
        self._nrps = defaultdict(dict)  # name -> method -> Prediction  # type: Dict[str, Dict[str, Prediction]]
        self.consensus = {}  # type: Dict[str, str]
        self.cluster_predictions = {}  # type: Dict[int, Tuple[str, bool]]
        self.domain_predictions = {}  # name -> Predictions # type: Dict[str, List[Predictions]]
        self.consensus_transat = {}  # type: Dict[str, str]

    @property
    def pks(self) -> PKSResults:
        """ The result of the PKS analyses """
        return self._pks

    @pks.setter
    def pks(self, value: PKSResults) -> None:
        assert isinstance(value, PKSResults)
        self._pks = value

    @property
    def nrps(self) -> Dict[str, Dict[str, Prediction]]:
        return self._nrps

    def to_json(self) -> Dict[str, Any]:
        results = {"schema_version": self._schema_version,
                   "record_id": self.record_id,
                   "pks": self.pks.to_json(),
                   "nrps": {},
                   "consensus": self.consensus,
                   "cluster_predictions": self.cluster_predictions}
        for domain, predictions in self.nrps.items():
            results["nrps"][domain] = {method: val.to_json() for method, val in predictions.items()},
        return results

    @staticmethod
    def from_json(json: Dict[str, Any], _record: Record) -> "NRPS_PKS_Results":
        assert "record_id" in json
        if json.get("schema_version") != NRPS_PKS_Results._schema_version:
            logging.warning("Mismatching schema version, dropping results")
            return None
        results = NRPS_PKS_Results(json["record_id"])
        pks = json.get("pks")
        if pks:
            results.pks = PKSResults.from_json(pks)
        nrps = json.get("nrps", {})
        for key, val in nrps.items():
            results.nrps[key] = PredictorSVMResult.from_json(val)  # TODO: unbreak
        results.consensus = json["consensus"]
        results.cluster_predictions.update(json["cluster_predictions"])
        return results

    def add_to_record(self, record: Record) -> None:
        """ Save substrate specificity predictions in NRPS/PKS domain sec_met info of record
        """

        for cds_feature in record.get_nrps_pks_cds_features():
            nrps_qualifier = cds_feature.nrps_pks
            for domain in nrps_qualifier.domains:
                feature = record.get_domain_by_name(domain.feature_name)
                assert isinstance(feature, AntismashDomain)
                domain.predictions.clear()
                if domain.name == "AMP-binding":
                    nrps_pred = generate_nrps_consensus(self.nrps[domain.feature_name])
                    domain.predictions["consensus"] = nrps_pred
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
