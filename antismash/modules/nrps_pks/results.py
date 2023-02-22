# License: GNU Affero General Public License v3 or later
# A copy of GNU AGPL v3 should have been included in this software package in LICENSE.txt.

""" Contains the results classes for the nrps_pks module """

from collections import defaultdict
import logging
from typing import Any, Dict, List, Optional

from antismash.common.module_results import ModuleResults
from antismash.common.secmet import Module, Record
from antismash.common.secmet.qualifiers import NRPSPKSQualifier
from antismash.detection.nrps_pks_domains import ModularDomain

from .parsers import generate_nrps_consensus
from .data_structures import Prediction, SimplePrediction
from .minowa.base import MinowaPrediction
from .nrpys import PredictorSVMResult
from .pks_names import get_short_form
from .at_analysis.at_analysis import ATPrediction


UNKNOWN = "(unknown)"


def modify_substrate(module: Module, base: str = "") -> str:  # pylint: disable=too-many-branches
    """ Builds a monomer including modifications from the given base.

        Arguments:
            module: the Module holding the relevant domains
            base: a string of the substate (or an empty string in case of trans-AT)

        Returns:
            the modified substrate, or an empty string if no appropriate base was
            given
    """
    if not module.is_complete():
        return ""

    domains = {domain.domain for domain in module.domains}

    if "PKS_KS" in domains and "PKS_AT" not in domains:
        base = "mal"

    if not base:
        return ""

    if module.module_type == Module.types.PKS:
        if "PKS_KR" in domains:
            conversions = {"mal": "ohmal", "mmal": "ohmmal", "mxmal": "ohmxmal", "emal": "ohemal"}
            base = conversions.get(base, base)

        if {"PKS_DH", "PKS_DH2", "PKS_DHt"}.intersection(domains):
            conversions = {"ohmal": "ccmal", "ohmmal": "ccmmal", "ohmxmal": "ccmxmal", "ohemal": "ccemal"}
            base = conversions.get(base, base)

        if "PKS_ER" in domains:
            conversions = {"ccmal": "redmal", "ccmmal": "redmmal", "ccmxmal": "redmxmal", "ccemal": "redemal"}
            base = conversions.get(base, base)

    state: List[str] = []
    for domain in module.domains:
        assert isinstance(domain, ModularDomain)
        if domain.domain != "MT":
            continue
        if domain.domain_subtype == "nMT":
            state.append("NMe")
        elif domain.domain_subtype == "cMT":
            state.append("Me")
        elif domain.domain_subtype == "oMT":
            state.append("OMe")

    if base.endswith("mmal"):
        state.append("Me")
        base = base.replace("mmal", "mal", 1)

    state.append(base)

    if "Epimerization" in domains:
        state.insert(0, "D")
    return "-".join(state)


class CandidateClusterPrediction:
    """ Contains prediction information for a CandidateCluster """
    def __init__(self, candidate_cluster_number: int, polymer: str,
                 domain_docking_used: bool, smiles: str = "", ordering: List[str] = None) -> None:
        self.candidate_cluster_number = candidate_cluster_number
        self.polymer = polymer
        self.domain_docking_used = domain_docking_used
        self.smiles = smiles
        self.ordering = ordering

    def to_json(self) -> Dict[str, Any]:
        """ Creates a JSON representation of a CandidateClusterPrediction """
        return {
            "sc_number": self.candidate_cluster_number,
            "polymer": self.polymer,
            "docking_used": self.domain_docking_used,
            "smiles": self.smiles,
            "ordering": self.ordering,
        }

    @staticmethod
    def from_json(json: Dict[str, Any]) -> "CandidateClusterPrediction":
        """ Rebuilds a CandidateClusterPrediction from a JSON dictionary """
        return CandidateClusterPrediction(json["sc_number"], json["polymer"],
                                          json["docking_used"], json["smiles"],
                                          json.get("ordering"))


class NRPS_PKS_Results(ModuleResults):
    """ The combined results of the nrps_pks module """
    schema_version = 3
    __slots__ = ["consensus", "consensus_transat", "region_predictions", "domain_predictions"]

    def __init__(self, record_id: str) -> None:
        super().__init__(record_id)
        # keep a mapping of domain name -> method -> Prediction
        self.domain_predictions: Dict[str, Dict[str, Prediction]] = defaultdict(dict)
        self.consensus: Dict[str, str] = {}  # domain name -> consensus
        self.region_predictions: Dict[int, List[CandidateClusterPrediction]] = defaultdict(list)
        self.consensus_transat: Dict[str, str] = {}

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
        domain_predictions: Dict[str, Dict[str, Any]] = defaultdict(dict)
        for domain, predictions in self.domain_predictions.items():
            domain_predictions[domain] = {method: val.to_json() for method, val in predictions.items()}
        region_json = {}
        for region, preds in self.region_predictions.items():
            region_json[region] = [pred.to_json() for pred in preds]
        results = {"schema_version": self.schema_version,
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
        if json.get("schema_version") != NRPS_PKS_Results.schema_version:
            logging.warning("Mismatching schema version, dropping results")
            return None
        results = NRPS_PKS_Results(json["record_id"])
        predictions = json.get("domain_predictions", {})
        for domain_name, method_predictions in predictions.items():
            for method, prediction in method_predictions.items():
                if method == "nrpys":
                    rebuilt: Prediction = PredictorSVMResult.from_json(prediction)
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
                prediction = CandidateClusterPrediction.from_json(pred)
                # the int conversion is important, since json can't use int keys
                results.region_predictions[int(region_number)].append(prediction)
        return results

    def _annotate_a_domain(self, domain: NRPSPKSQualifier.Domain) -> None:
        assert domain.name in ["AMP-binding", "A-OX"]
        predictions = self.domain_predictions[domain.feature_name]
        domain.add_prediction("substrate consensus", generate_nrps_consensus(predictions))

    def _annotate_at_domain(self, domain: NRPSPKSQualifier.Domain, transat_cluster: bool) -> None:
        assert domain.name == "PKS_AT"
        predictions = self.domain_predictions[domain.feature_name]

        if transat_cluster:
            consensus = self.consensus_transat[domain.feature_name]
        else:
            consensus = self.consensus[domain.feature_name]
        domain.add_prediction("substrate consensus", consensus)

        sig = UNKNOWN
        pks_sig = predictions["signature"]
        if pks_sig and pks_sig.get_classification():
            sig = pks_sig.get_classification()[0]
        domain.add_prediction("PKS signature", sig)

        minowa = predictions["minowa_at"].get_classification()[0]
        domain.add_prediction("Minowa", get_short_form(minowa))

    def _annotate_cal_domain(self, domain: NRPSPKSQualifier.Domain) -> None:
        assert domain.name == "CAL_domain"
        minowa = self.domain_predictions[domain.feature_name]["minowa_cal"].get_classification()[0]
        domain.add_prediction("Minowa", get_short_form(minowa))

    def _annotate_kr_domain(self, domain: NRPSPKSQualifier.Domain) -> None:
        assert domain.name == "PKS_KR"
        predictions = self.domain_predictions[domain.feature_name]

        activity = predictions["kr_activity"].get_classification()[0]
        domain.add_prediction("KR activity", activity)

        stereo_chem = UNKNOWN
        stereo_pred = predictions.get("kr_stereochem")
        if stereo_pred:
            stereo_chem = stereo_pred.get_classification()[0]
        domain.add_prediction("KR stereochemistry", stereo_chem)

    def add_to_record(self, record: Record) -> None:
        """ Save substrate specificity predictions in NRPS/PKS domain sec_met info of record
        """
        for candidate_cluster_preds in self.region_predictions.values():
            for cluster_pred in candidate_cluster_preds:
                assert isinstance(cluster_pred, CandidateClusterPrediction), type(cluster_pred)
                candidate = record.get_candidate_cluster(cluster_pred.candidate_cluster_number)
                candidate.smiles_structure = cluster_pred.smiles

        for cds_feature in record.get_nrps_pks_cds_features():
            assert cds_feature.region, "CDS parent region removed since analysis"
            nrps_qualifier = cds_feature.nrps_pks
            for domain in nrps_qualifier.domains:
                feature = record.get_domain_by_name(domain.feature_name)
                assert isinstance(feature, ModularDomain)

                if domain.name in ["AMP-binding", "A-OX"]:
                    self._annotate_a_domain(domain)
                elif domain.name == "PKS_AT":
                    self._annotate_at_domain(domain, "transatpks" in cds_feature.region.products)
                elif domain.name == "CAL_domain":
                    self._annotate_cal_domain(domain)
                elif domain.name == "PKS_KR":
                    self._annotate_kr_domain(domain)
                # otherwise one of many without prediction methods/relevance (PCP, Cglyc, etc)

                for method, pred in domain.get_predictions().items():
                    feature.specificity.append(f"{method}: {pred}")

            for module in cds_feature.modules:
                if not module.is_complete():
                    continue
                substrate = ""
                for module_domain in module.domains:
                    consensus = self.consensus.get(module_domain.get_name())
                    if consensus:
                        substrate = consensus
                        break
                if not substrate:  # probably a trans-AT PKS
                    domains = {dom.domain for dom in module.domains}
                    if module.module_type == module.types.PKS and "PKS_AT" not in domains:
                        substrate = "mal"
                    else:
                        raise ValueError(f"missing substrate in non-transAT module: {module}")

                module.add_monomer(substrate, modify_substrate(module, substrate))
