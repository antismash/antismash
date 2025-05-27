# License: GNU Affero General Public License v3 or later
# A copy of GNU AGPL v3 should have been included in this software package in LICENSE.txt.

""" Contains the results classes for the terpene module """

import logging
from typing import Any, Optional, Self
from dataclasses import asdict, dataclass

from antismash.common.module_results import ModuleResults
from antismash.common.secmet import Record

from .data_loader import (
    CompoundGroup,
    Reaction,
    load_compounds,
    MissingCompoundError,
    MissingHmmError,
)


@dataclass(slots=True, eq=True)
class DomainPrediction:
    """ A prediction for a terpene biosynthetic domain
    """
    domain_type: str
    subtypes: tuple[str, ...]
    start: int
    end: int
    reactions: tuple[Reaction, ...]

    def __str__(self) -> str:
        return (
            f"DomainPrediction(type={self.domain_type}, subtypes={self.subtypes}, start={self.start}, "
            f"end={self.end}, reactions={[str(reaction) for reaction in self.reactions]})"
        )

    def to_json(self) -> dict[str, Any]:
        """ Returns a JSON-friendly representation """
        data = asdict(self)
        reactions = data.pop("reactions")
        if reactions:
            data["reactions"] = [reaction.to_json() for reaction in reactions]
        return data

    @classmethod
    def from_json(cls, data: dict[str, Any], compound_groups: dict[str, CompoundGroup],
                  ) -> Self:
        """ Reconstructs an instance from a JSON representation """
        reactions = tuple(Reaction.from_json(reaction, compound_groups) for reaction in data.pop("reactions", []))
        subtypes = tuple(data.pop("subtypes"))
        return cls(**data, subtypes=subtypes, reactions=reactions)


class ProtoclusterPrediction:
    """ A prediction for a terpene protocluster
    """
    def __init__(self, cds_predictions: dict[str, list[DomainPrediction]],
                 products: list[CompoundGroup] = None) -> None:
        self.cds_predictions = cds_predictions
        self.products = products or []

    def add_product(self, product: CompoundGroup) -> None:
        """ Adds a product to products """
        self.products.append(product)

    def get_unique_compounds(self) -> set[CompoundGroup]:
        """ Retrieve all compounds associated with this ProtoclusterPrediction """
        compounds: set[CompoundGroup] = set()
        for domains in self.cds_predictions.values():
            for domain in domains:
                for reaction in domain.reactions:
                    compounds.update(reaction.substrates + reaction.products)
        return compounds

    def get_functional_groups(self) -> set[str]:
        """ Retrieve functional groups of the products """
        func_groups: set[str] = set()
        for product in self.products:
            func_groups.update(product.functional_groups)
        return func_groups

    def __eq__(self, other: Any) -> bool:
        if not isinstance(other, ProtoclusterPrediction):
            return False
        return self.cds_predictions == other.cds_predictions and self.products == other.products

    def __str__(self) -> str:
        parts = [
            f"CDSs: {len(self.cds_predictions)}",
        ]
        for cds, predictions in self.cds_predictions.items():
            parts.append(str(cds))
            parts.append(("\n ".join(map(str, predictions))))
        return "\n".join(parts)

    def to_json(self) -> dict[str, Any]:
        """ Returns a JSON-friendly representation """
        cds_preds = {}
        for name, preds in self.cds_predictions.items():
            cds_preds[name] = [pred.to_json() for pred in preds]
        return {
            "cds_predictions": cds_preds,
            "products": [product.name for product in self.products],
        }

    @classmethod
    def from_json(cls, data: dict[str, Any], compound_groups: dict[str, CompoundGroup]
                  ) -> Self:
        """ Reconstructs an instance from a JSON representation """
        cds_preds: dict[str, list[DomainPrediction]] = {}
        for name, preds in data.pop("cds_predictions").items():
            domains = [DomainPrediction.from_json(pred, compound_groups) for pred in preds]
            cds_preds[name] = domains
        try:
            return cls(cds_preds, [compound_groups[compound] for compound in data["products"]])
        except KeyError as key:
            raise MissingCompoundError(f"Compound group {key} is not defined.")


class TerpeneResults(ModuleResults):
    """ The combined results of the terpene module """
    _schema_version = 1
    __slots__ = ["cluster_predictions"]

    def __init__(self, record_id: str) -> None:
        super().__init__(record_id)
        self.cluster_predictions: dict[int, ProtoclusterPrediction] = {}

    def __repr__(self) -> str:
        return f"TerpeneResults(clusters={list(self.cluster_predictions)})"

    def __str__(self) -> str:
        parts = []
        for cluster_id, prediction in self.cluster_predictions.items():
            parts.append(f"Protocluster {cluster_id}\n")
            parts.append(str(prediction))
        return "".join(parts)

    def to_json(self) -> dict[str, Any]:
        """ Returns a JSON-friendly representation """
        predictions = {cluster_number: pred.to_json() for cluster_number, pred in self.cluster_predictions.items()}
        results = {
            "schema_version": self._schema_version,
            "record_id": self.record_id,
            "protocluster_predictions": predictions,
        }
        return results

    @staticmethod
    def from_json(data: dict[str, Any], record: Record,
                  ) -> Optional["TerpeneResults"]:
        """ Reconstructs an instance from a JSON representation """
        assert "record_id" in data
        if data.get("schema_version") != TerpeneResults._schema_version:
            logging.warning("Mismatching schema version, dropping terpene results")
            return None

        results = TerpeneResults(data["record_id"])
        try:
            for cluster_id, prediction in data["protocluster_predictions"].items():
                cluster_prediction = ProtoclusterPrediction.from_json(prediction, load_compounds())
                results.cluster_predictions[int(cluster_id)] = cluster_prediction
        except (MissingCompoundError, MissingHmmError) as err:
            logging.warning("Discarding results, missing referenced terpene data: %s", err)
            return None
        return results

    def add_to_record(self, record: Record) -> None:
        # nothing to store
        pass
