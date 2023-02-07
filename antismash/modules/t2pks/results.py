# License: GNU Affero General Public License v3 or later
# A copy of GNU AGPL v3 should have been included in this software package in LICENSE.txt.

""" Contains the results classes for the nrps_pks module """

import logging
from typing import Any, Dict, List, Optional, Union
from typing import Set, Tuple

from antismash.common.module_results import ModuleResults
from antismash.common.secmet import Record, GeneFunction
from antismash.common.secmet.qualifiers import T2PKSQualifier


class Prediction:
    """ A generic named score and evalue pair that allows JSON conversion
    """
    def __init__(self, name: str, score: float, evalue: float) -> None:
        self.name = str(name)
        self.score = float(score)
        self.evalue = float(evalue)

    def __repr__(self) -> str:
        return f"Prediction({self.name}, {self.score:.1f}, {self.evalue:g})"

    def __str__(self) -> str:
        return f"{self.name} (Score: {self.score:.1f}; E-value: {self.evalue:g})"

    def __eq__(self, other: Any) -> bool:
        if not isinstance(other, Prediction):
            return False
        return (self.name == other.name
                and self.score == other.score
                and self.evalue == other.evalue)

    def to_json(self) -> Tuple[str, float, float]:
        """ Returns a JSON-friendly representation of the Prediction """
        return self.name, self.score, self.evalue

    @staticmethod
    def from_json(json: List[Union[str, float]]) -> "Prediction":
        """ Reconstructs a Prediction from a JSON representation """
        return Prediction(str(json[0]), float(json[1]), float(json[2]))


class CDSPrediction:
    """ Contains the CDS-specific predictions for the type II PKS module """
    def __init__(self, protein_type: str, protein_function: Optional[str],
                 bitscore: float, evalue: float) -> None:
        self.ptype = protein_type
        self.pfunc = protein_function
        self.bitscore = bitscore
        self.evalue = evalue

    def __repr__(self) -> str:
        return f"CDSPrediction({self.ptype}, {self.pfunc})"

    def __str__(self) -> str:
        base = self.ptype
        if self.pfunc is not None:
            base = f"{self.ptype} {self.pfunc}"
        return f"{base} (Score: {self.bitscore:.1f}; E-value: {self.evalue:g})"

    def to_json(self) -> Tuple[str, Optional[str], float, float]:
        """ Converts a CDSPrediction into a JSON friendly format """
        return (self.ptype, self.pfunc, self.bitscore, self.evalue)

    @staticmethod
    def from_json(json: Tuple[str, Optional[str], float, float]) -> "CDSPrediction":
        """ Rebuilds a CDSPrediction from JSON """
        if len(json) != 4:
            raise ValueError(f"Invalid CDSPrediction JSON: {json}")
        func = None
        if json[1] is not None:
            func = str(json[1])
        return CDSPrediction(str(json[0]), func, float(json[2]), float(json[3]))


class ProtoclusterPrediction:
    """ A prediction for a single Cluster, including starter units, elongations,
        weights and classes.
    """
    def __init__(self, cds_predictions: Dict[str, List[CDSPrediction]],  # pylint: disable=too-many-arguments
                 starter_units: List[Prediction],
                 malonyl_elongations: List[Prediction],
                 product_classes: Set[str],
                 molecular_weights: Dict[str, float],
                 start: int,
                 end: int) -> None:
        self.cds_predictions = cds_predictions
        self.starter_units = starter_units
        self.malonyl_elongations = malonyl_elongations
        self.product_classes = product_classes
        self.molecular_weights = molecular_weights
        self.start = start
        self.end = end

    def __repr__(self) -> str:
        starters = [unit.name for unit in self.starter_units]
        elongations = [elong.name for elong in self.malonyl_elongations]
        classes = sorted(self.product_classes)
        return f"ClusterPrediction(starters={starters}, elongations={elongations}, classes={classes})"

    def __str__(self) -> str:
        parts = [
            f"Starter units: {self.starter_units}",
            f"Malonyl elongations: {self.malonyl_elongations}",
            f"Product classes: {','.join(sorted(self.product_classes))}",
            f"Molecular weights: {self.molecular_weights}",
            f"CDSs: {len(self.cds_predictions)}",
        ]

        for cds, predictions in self.cds_predictions.items():
            parts.append(str(cds))
            parts.append(" " + ("\n ".join(map(str, predictions))))
        return "\n".join(parts)

    def to_json(self) -> Dict[str, Any]:
        """ Converts a ClusterPrediction into a JSON friendly format """
        cds_preds = {}
        for name, preds in self.cds_predictions.items():
            cds_preds[name] = [pred.to_json() for pred in preds]

        return {"starter_units": [pred.to_json() for pred in self.starter_units],
                "elongations": [pred.to_json() for pred in self.malonyl_elongations],
                "product_classes": list(self.product_classes),
                "mol_weights": self.molecular_weights,
                "cds_preds": cds_preds,
                "start": self.start,
                "end": self.end
                }

    @staticmethod
    def from_json(json: Dict[str, Any]) -> "ProtoclusterPrediction":
        """ Rebuilds a ProtoclusterPrediction from JSON """
        assert isinstance(json, dict), json
        cds_predictions = {name: list(map(CDSPrediction.from_json, preds)) for name, preds in json["cds_preds"].items()}
        # JSON doesn't do tuples, so convert back to tuples
        starters = [Prediction.from_json(start) for start in json["starter_units"]]
        elongations = [Prediction.from_json(elong) for elong in json["elongations"]]
        return ProtoclusterPrediction(cds_predictions, starters, elongations,
                                      set(json["product_classes"]), json["mol_weights"],
                                      int(json["start"]), int(json["end"]))


class T2PKSResults(ModuleResults):
    """ The combined results of the type 2 PKS module """
    _schema_version = 3
    __slots__ = ["cluster_predictions"]

    def __init__(self, record_id: str) -> None:
        super().__init__(record_id)
        self.cluster_predictions: Dict[int, ProtoclusterPrediction] = {}

    def __repr__(self) -> str:
        return f"T2PKSResults(clusters={list(self.cluster_predictions)})"

    def __str__(self) -> str:
        parts = []
        for cluster_id, prediction in self.cluster_predictions.items():
            parts.append(f"Protocluster {cluster_id}\n")
            parts.append(str(prediction))
        return "".join(parts)

    def to_json(self) -> Dict[str, Any]:
        clusters = {cluster_number: pred.to_json() for cluster_number, pred in self.cluster_predictions.items()}
        results = {"schema_version": self._schema_version,
                   "record_id": self.record_id,
                   "protocluster_predictions": clusters}
        return results

    @staticmethod
    def from_json(json: Dict[str, Any], _record: Record) -> Optional["T2PKSResults"]:
        assert "record_id" in json
        if json.get("schema_version") != T2PKSResults._schema_version:
            logging.warning("Mismatching schema version, dropping T2PKS results")
            return None
        results = T2PKSResults(json["record_id"])

        for cluster_id, json_prediction in json["protocluster_predictions"].items():
            cluster_prediction = ProtoclusterPrediction.from_json(json_prediction)
            # int is required because JSON keys are strings
            results.cluster_predictions[int(cluster_id)] = cluster_prediction

        return results

    def add_to_record(self, record: Record) -> None:
        """ Save type II PKS prediction in record.

            Cluster predictions are saved the relevant Cluster feature.
            Gene functions are added to each CDSFeature.
        """
        for cluster_number, prediction in self.cluster_predictions.items():
            cluster = record.get_protocluster(cluster_number)
            cluster.t2pks = T2PKSQualifier(list(map(str, prediction.starter_units)),
                                           list(map(str, prediction.malonyl_elongations)),
                                           sorted(prediction.product_classes),
                                           prediction.molecular_weights)

            for cds, sub_predictions in prediction.cds_predictions.items():
                for pred in sub_predictions:
                    cds_feature = record.get_cds_by_name(cds)
                    cds_feature.gene_functions.add(GeneFunction.ADDITIONAL, 't2pks', str(pred))
