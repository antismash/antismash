# License: GNU Affero General Public License v3 or later
# A copy of GNU AGPL v3 should have been included in this software package in LICENSE.txt.

""" A class for protocluster features """

from typing import Any, Dict, List, Optional, Set, Type, TypeVar

from Bio.SeqFeature import SeqFeature

from .cds_feature import CDSFeature
from .cdscollection import CDSCollection, CoredCollectionMixin
from .feature import Feature, FeatureLocation
from ..locations import Location, location_from_string
from ..qualifiers.t2pks import T2PKSQualifier
from ..qualifiers.gene_functions import GeneFunction

T = TypeVar("T", bound="Protocluster")
S = TypeVar("S", bound="SideloadedProtocluster")


class Protocluster(CDSCollection, CoredCollectionMixin):
    """ A feature which marks a specific region of a record as interesting.
        Protoclusters are only those determined by a rule-based method of detection
        and with a defined product.

        The core location covers the first and last CDS feature that caused the
        protocluster to be formed, while the surrounding location includes the context.
    """
    core_seqfeature_type = "proto_core"
    __slots__ = ["_core_location", "detection_rule", "product", "product_category",
                 "tool", "cutoff",
                 "_definition_cdses", "neighbourhood_range", "t2pks",]
    FEATURE_TYPE = "protocluster"  # primary type only

    def __init__(self, core_location: Location, surrounding_location: Location,
                 tool: str, product: str, cutoff: int, neighbourhood_range: int,
                 detection_rule: str, product_category: str = "other") -> None:

        if core_location.crosses_origin() and not surrounding_location.crosses_origin():
            raise ValueError(f"a core location ({core_location}) crossing the origin requires "
                             f"the surrounding area ({surrounding_location}) to also cross the origin")
        assert len(surrounding_location.parts) >= len(core_location.parts)
        super().__init__(surrounding_location, feature_type=self.FEATURE_TYPE)
        # cluster-wide
        self.detection_rule = detection_rule
        if not product.replace("-", "").replace("_", "").isalnum() or product[0] in "-_" or product[-1] in "-_":
            raise ValueError(f"invalid protocluster product: {product}")
        self.product = product
        self.product_category = product_category
        self.tool = tool

        # core specific
        self._core_location = core_location
        self.cutoff = cutoff
        self._definition_cdses: Set[CDSFeature] = set()

        # neighbourhood specific
        self.neighbourhood_range = neighbourhood_range

        # analysis annotations
        self.t2pks: Optional[T2PKSQualifier] = None

    def __str__(self) -> str:
        return f"Protocluster({self.location}, product={self.product})"

    def get_protocluster_number(self) -> int:
        """ Returns the protoclusters's numeric ID, only guaranteed to be consistent
            when the same protoclusters are defined in the parent record
        """
        if not self._parent_record:
            raise ValueError("Protocluster not in a record")
        return self._parent_record.get_protocluster_number(self)

    @property
    def contig_edge(self) -> bool:
        # always trust an explicit positive
        if super().contig_edge:
            return True
        # then check, and update the stored value, if either cutoff or neighbourhood extend that far
        start = min(self.location.start, self.core_location.start - self.cutoff)
        end = max(self.location.end, self.core_location.end + self.cutoff)
        contig_edge = start < 0 or end >= len(self.parent_record.seq)
        self._contig_edge = contig_edge
        return contig_edge

    @property
    def core_location(self) -> Location:
        return self._core_location

    @property
    def definition_cdses(self) -> Set[CDSFeature]:
        """ Returns the set of CDSFeatures responsible for the creation of this protocluster """
        return set(self._definition_cdses)

    def add_cds(self, cds: CDSFeature) -> None:
        super().add_cds(cds)
        if not cds.is_contained_by(self.core_location):
            return
        cores = cds.gene_functions.get_by_function(GeneFunction.CORE)
        if any(core.product == self.product for core in cores):
            self._definition_cdses.add(cds)

    def core_crosses_origin(self) -> bool:
        """ Returns True if the core of the protocluster crosses the origin.
        """
        return len(self.core_location.parts) > 1

    def to_biopython(self, qualifiers: Optional[Dict[str, List[str]]] = None) -> List[SeqFeature]:
        common = {
            "neighbourhood": [str(self.neighbourhood_range)],
            "cutoff": [str(self.cutoff)],
            "product": [self.product],
            "aStool": [self.tool],
            "detection_rule": [self.detection_rule]
        }
        if self._parent_record:
            common["protocluster_number"] = [str(self.get_protocluster_number())]

        shared_qualifiers = dict(qualifiers) if qualifiers else {}
        shared_qualifiers.update(common)
        if self.t2pks:
            shared_qualifiers.update(self.t2pks.to_biopython_qualifiers())

        core_qualifiers = dict(self._qualifiers)
        core_qualifiers.update(shared_qualifiers)
        core_feature = SeqFeature(self.core_location, type=self.core_seqfeature_type)
        for key, val in sorted(core_qualifiers.items()):
            core_feature.qualifiers[key] = val
            # this doesn't go through to the Feature annotations, so add the tool explicitly
            core_feature.qualifiers["tool"] = ["antismash"]

        shared_qualifiers["core_location"] = [str(self.core_location)]
        if self.product_category:
            shared_qualifiers["category"] = [self.product_category]
        neighbourhood_feature = super().to_biopython(shared_qualifiers)[0]
        return [neighbourhood_feature, core_feature]

    @classmethod
    def from_biopython(cls: Type[T], bio_feature: SeqFeature, feature: T = None,
                       leftovers: Dict[str, List[str]] = None, record: Any = None) -> T:
        assert bio_feature.type == Protocluster.FEATURE_TYPE
        if leftovers is None:
            leftovers = Feature.make_qualifiers_copy(bio_feature)
        if leftovers.get("aStool", [""])[0].startswith("externally annotated"):
            external = SideloadedProtocluster.from_biopython(bio_feature)
            assert isinstance(external, cls)
            return external

        if not feature:
            category = leftovers.pop("category", [""])[0]
            # grab mandatory qualifiers and create the class
            try:
                neighbourhood_range = int(leftovers.pop("neighbourhood")[0])
                cutoff = int(leftovers.pop("cutoff")[0])
                product = leftovers.pop("product")[0]
                tool = leftovers.pop("aStool")[0]
                rule = leftovers.pop("detection_rule")[0]
                core_location = location_from_string(leftovers.pop("core_location")[0])
            except KeyError as err:
                raise ValueError(f"{cls.FEATURE_TYPE} missing expected qualifier: {err}")
            feature = cls(core_location, bio_feature.location,
                          tool, product, cutoff, neighbourhood_range, rule, product_category=category)

        # remove run-specific info
        leftovers.pop("protocluster_number", "")

        # rebuild analysis annotations
        feature.t2pks = T2PKSQualifier.from_biopython_qualifiers(leftovers)

        # grab optional parent qualifiers
        updated = super().from_biopython(bio_feature, feature, leftovers, record=record)
        assert updated is feature, f"feature changed: {feature} -> {updated}"
        assert isinstance(updated, Protocluster)
        return updated


class SideloadedProtocluster(Protocluster):
    """ A variant of Protocluster specifically for sideloaded features
    """
    __slots__ = ["extra_qualifiers"]

    def __init__(self, core_location: FeatureLocation, surrounding_location: FeatureLocation,
                 tool: str, product: str, neighbourhood_range: int = 0,
                 extra_qualifiers: Dict[str, List[str]] = None) -> None:
        if not neighbourhood_range:
            neighbourhood_range = max(core_location.start - surrounding_location.start,
                                      surrounding_location.end - core_location.end)
        super().__init__(core_location, surrounding_location, tool, product,
                         cutoff=0, neighbourhood_range=neighbourhood_range,
                         detection_rule="from external annotation")
        self.extra_qualifiers = extra_qualifiers or {}

    @property
    def definition_cdses(self) -> Set[CDSFeature]:
        # a sideloaded protocluster cannot have definition cdses
        return set()

    def to_biopython(self, qualifiers: Optional[Dict[str, List[str]]] = None) -> List[SeqFeature]:
        features = super().to_biopython(qualifiers)
        for feature in features:
            feature.qualifiers["aStool"] = [f"externally annotated by: {self.tool}"]
            if self.extra_qualifiers:
                feature.qualifiers["external_qualifier_ids"] = list(self.extra_qualifiers)
                feature.qualifiers.update(self.extra_qualifiers)
        return features

    @classmethod
    def from_biopython(cls: Type[S], bio_feature: SeqFeature, feature: S = None,
                       leftovers: Optional[Dict] = None, record: Any = None) -> S:
        if leftovers is None:
            leftovers = Feature.make_qualifiers_copy(bio_feature)
        if not feature:
            tool = leftovers.pop("aStool")[0]
            tool = tool.split(": ", 1)[1]
            leftovers["aStool"] = [tool]
            core_location = location_from_string(leftovers["core_location"][0])
            product = leftovers.pop("product")[0]
            extra_qualifiers: Dict[str, List[str]] = {}
            for key in leftovers.pop("external_qualifier_ids", []):
                if key in leftovers:
                    extra_qualifiers[key] = leftovers.pop(key)
            neighbourhood_range = int(leftovers.pop("neighbourhood")[0])
            feature = cls(core_location, bio_feature.location, tool,
                          product, neighbourhood_range=neighbourhood_range,
                          extra_qualifiers=extra_qualifiers)
        assert isinstance(feature, cls)
        feature = super().from_biopython(bio_feature, feature=feature, leftovers=leftovers)
        return feature
