# License: GNU Affero General Public License v3 or later
# A copy of GNU AGPL v3 should have been included in this software package in LICENSE.txt.

""" Contains the results classes for the nrps_pks module """

from collections import defaultdict
import logging
from typing import Any, Dict, List, Optional  # pylint: disable=unused-import

from antismash.common.module_results import ModuleResults
from antismash.common.secmet import Record, AntismashDomain
from antismash.common.secmet.qualifiers import NRPSPKSQualifier

from .parsers import generate_nrps_consensus
from .data_structures import Prediction, SimplePrediction
from .minowa.base import MinowaPrediction
from .nrps_predictor import PredictorSVMResult
from .pks_names import get_short_form
from .at_analysis.at_analysis import ATPrediction

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


UNKNOWN = "(unknown)"


class SuperClusterPrediction:
    """ Contains prediction information for a SuperCluster """
    def __init__(self, supercluster_number: int, polymer: str,
                 domain_docking_used: bool, smiles: str = "") -> None:
        self.supercluster_number = supercluster_number
        self.polymer = polymer
        self.domain_docking_used = domain_docking_used
        self.smiles = smiles

    def to_json(self) -> Dict[str, Any]:
        """ Creates a JSON representation of a SuperClusterPrediction """
        return {
            "sc_number": self.supercluster_number,
            "polymer": self.polymer,
            "docking_used": self.domain_docking_used,
            "smiles": self.smiles
        }

    @staticmethod
    def from_json(json: Dict[str, Any]) -> "SuperClusterPrediction":
        """ Rebuilds a SuperClusterPrediction from a JSON dictionary """
        return SuperClusterPrediction(json["sc_number"], json["polymer"],
                                      json["docking_used"], json["smiles"])


class NRPS_PKS_Results(ModuleResults):
    """ The combined results of the nrps_pks module """
    _schema_version = 3
    __slots__ = ["consensus", "consensus_transat", "region_predictions", "domain_predictions"]

    def __init__(self, record_id: str) -> None:
        super().__init__(record_id)
        # keep a mapping of domain name -> method -> Prediction
        self.domain_predictions = defaultdict(dict)  # type: Dict[str, Dict[str, Prediction]]
        self.consensus = {}  # type: Dict[str, str]  # domain name -> consensus
        self.region_predictions = defaultdict(list)  # type: Dict[int, List[SuperClusterPrediction]]
        self.consensus_transat = {}  # type: Dict[str, str]

    def add_method_results(self, method: str, results: Dict[str, Prediction]) -> None:
        """ Add per-domain results for a single prediction method

            Arguments:
                method: the name of the method that generated the results
                results: a dictionary mapping the domain name to a Prediction subclass

            Returns:
                None
        """
        for domain_name, prediction in results.items():
            self.domain_predictions[domain_name][method] = prediction

    def to_json(self) -> Dict[str, Any]:
        domain_predictions = defaultdict(dict)  # type: Dict[str, Dict[str, Any]]
        for domain, predictions in self.domain_predictions.items():
            domain_predictions[domain] = {method: val.to_json() for method, val in predictions.items()}
        region_json = {}
        for region, preds in self.region_predictions.items():
            region_json[region] = [pred.to_json() for pred in preds]
        results = {"schema_version": self._schema_version,
                   "record_id": self.record_id,
                   "domain_predictions": domain_predictions,
                   "consensus": self.consensus,
                   "consensus_transat": self.consensus_transat,
                   "region_predictions": region_json,
                   }
        return results

    @staticmethod
    def from_json(json: Dict[str, Any], _record: Record) -> Optional["NRPS_PKS_Results"]:
        assert "record_id" in json
        if json.get("schema_version") != NRPS_PKS_Results._schema_version:
            logging.warning("Mismatching schema version, dropping results")
            return None
        results = NRPS_PKS_Results(json["record_id"])
        predictions = json.get("domain_predictions", {})
        for domain_name, method_predictions in predictions.items():
            for method, prediction in method_predictions.items():
                if method == "NRPSPredictor2":
                    rebuilt = PredictorSVMResult.from_json(prediction)  # type: Prediction
                elif method.startswith("minowa"):
                    rebuilt = MinowaPrediction.from_json(prediction)
                elif method == "signature":
                    rebuilt = ATPrediction.from_json(prediction)
                else:
                    rebuilt = SimplePrediction.from_json(prediction)
                results.domain_predictions[domain_name][method] = rebuilt
        results.consensus = json["consensus"]
        results.consensus_transat = json["consensus_transat"]
        for region_number, predictions in json["region_predictions"].items():
            for pred in predictions:
                prediction = SuperClusterPrediction.from_json(pred)
                # the int conversion is important, since json can't use int keys
                results.region_predictions[int(region_number)].append(prediction)
        return results

    def _annotate_a_domain(self, domain: NRPSPKSQualifier.Domain) -> None:
        assert domain.name in ["AMP-binding", "A-OX"]
        predictions = self.domain_predictions[domain.feature_name]
        domain.predictions["consensus"] = generate_nrps_consensus(predictions)

    def _annotate_at_domain(self, domain: NRPSPKSQualifier.Domain, transat_cluster: bool) -> None:
        assert domain.name == "PKS_AT"
        predictions = self.domain_predictions[domain.feature_name]

        if transat_cluster:
            consensus = self.consensus_transat[domain.feature_name]
        else:
            consensus = self.consensus[domain.feature_name]
        domain.predictions["consensus"] = consensus

        sig = UNKNOWN
        pks_sig = predictions["signature"]
        if pks_sig and pks_sig.get_classification():
            sig = pks_sig.get_classification()[0]
        domain.predictions["PKS signature"] = sig

        minowa = predictions["minowa_at"].get_classification()[0]
        domain.predictions["Minowa"] = get_short_form(minowa)

    def _annotate_cal_domain(self, domain: NRPSPKSQualifier.Domain) -> None:
        assert domain.name == "CAL_domain"
        minowa = self.domain_predictions[domain.feature_name]["minowa_cal"].get_classification()[0]
        domain.predictions["Minowa"] = get_short_form(minowa)

    def _annotate_kr_domain(self, domain: NRPSPKSQualifier.Domain) -> None:
        assert domain.name == "PKS_KR"
        predictions = self.domain_predictions[domain.feature_name]

        activity = predictions["kr_activity"].get_classification()[0]
        domain.predictions["KR activity"] = activity

        stereo_chem = UNKNOWN
        stereo_pred = predictions.get("kr_stereochem")
        if stereo_pred:
            stereo_chem = stereo_pred.get_classification()[0]
        domain.predictions["KR stereochemistry"] = stereo_chem

    def add_to_record(self, record: Record) -> None:
        """ Save substrate specificity predictions in NRPS/PKS domain sec_met info of record
        """
        for supercluster_preds in self.region_predictions.values():
            for cluster_pred in supercluster_preds:
                assert isinstance(cluster_pred, SuperClusterPrediction), type(cluster_pred)
                record.get_supercluster(cluster_pred.supercluster_number).smiles_structure = cluster_pred.smiles

        for cds_feature in record.get_nrps_pks_cds_features():
            assert cds_feature.region, "CDS parent region removed since analysis"
            nrps_qualifier = cds_feature.nrps_pks
            for domain in nrps_qualifier.domains:
                feature = record.get_domain_by_name(domain.feature_name)
                assert isinstance(feature, AntismashDomain)

                domain.predictions.clear()
                if domain.name in ["AMP-binding", "A-OX"]:
                    self._annotate_a_domain(domain)
                elif domain.name == "PKS_AT":
                    self._annotate_at_domain(domain, "transatpks" in cds_feature.region.products)
                elif domain.name == "CAL_domain":
                    self._annotate_cal_domain(domain)
                elif domain.name == "PKS_KR":
                    self._annotate_kr_domain(domain)
                # otherwise one of many without prediction methods/relevance (PCP, Cglyc, etc)

                for method, pred in domain.predictions.items():
                    feature.specificity.append("%s: %s" % (method, pred))

                mapping = DOMAIN_TYPE_MAPPING.get(domain.name)
                if mapping:
                    feature.domain_subtype = domain.name
                    feature.domain = mapping
