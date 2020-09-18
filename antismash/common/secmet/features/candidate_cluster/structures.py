# License: GNU Affero General Public License v3 or later
# A copy of GNU AGPL v3 should have been included in this software package in LICENSE.txt.

""" Contains the classes required for candidate clusters """

from collections import OrderedDict
from enum import Enum, unique
from typing import Any, Dict, List, Optional, Tuple, Type, TypeVar

from Bio.SeqFeature import SeqFeature

from ..cdscollection import CDSCollection
from ..protocluster import Protocluster
from ..feature import FeatureLocation, Feature
from ...locations import combine_locations

T = TypeVar("T", bound="CandidateCluster")


@unique
class CandidateClusterKind(Enum):
    """ An Enum representing the kind of a CandidateCluster.
        Allows for more flexible conversion and more robust value constraints.
    """
    SINGLE = 0
    INTERLEAVED = 1
    NEIGHBOURING = 2
    CHEMICAL_HYBRID = 3

    def __str__(self) -> str:
        return str(self.name).lower()

    @staticmethod
    def from_string(label: str) -> "CandidateClusterKind":
        """ Converts a string to a CandidateClusterKind instance when possible.
            Raises an error if not possible.
        """
        for value in CandidateClusterKind:
            if str(value) == label:
                return value
        raise ValueError("unknown candidate cluster kind: %s" % label)


class CandidateCluster(CDSCollection):
    """ A class representing a collection of overlapping Cluster features.
        The location of a CandidateCluster is defined as the minimum area that would
        contain all of the child Protolusters.
    """
    FEATURE_TYPE = "cand_cluster"
    kinds = CandidateClusterKind
    __slots__ = ["_protoclusters", "_kind", "_core_location", "smiles_structure", "polymer"]

    def __init__(self, kind: CandidateClusterKind, protoclusters: List[Protocluster],
                 smiles: str = None, polymer: str = None) -> None:
        if not protoclusters:
            raise ValueError("A CandidateCluster cannot exist without at least one Protocluster")
        for protocluster in protoclusters:
            assert isinstance(protocluster, Protocluster), type(protocluster)
        if not isinstance(kind, CandidateClusterKind):
            raise TypeError("argument 1 should be CandidateClusterKind, had %s" % type(kind))
        location = combine_locations(cluster.location for cluster in protoclusters)
        super().__init__(location, feature_type=CandidateCluster.FEATURE_TYPE, child_collections=protoclusters)
        self._protoclusters = protoclusters
        self._kind = kind
        self.smiles_structure = smiles
        self.polymer = polymer
        self._core_location = None

    def __repr__(self) -> str:
        return "CandidateCluster(%s, %s)" % (self.location, self.kind)

    def get_candidate_cluster_number(self) -> int:
        """ Returns the candidate clusters's numeric ID, only guaranteed to be consistent for
            when the same clusters and subregions are defined in the parent record
        """
        if not self._parent_record:
            raise ValueError("CandidateCluster not contained in record")
        return self._parent_record.get_candidate_cluster_number(self)

    @property
    def kind(self) -> CandidateClusterKind:
        """ The kind of CandidateCluster """
        return self._kind

    @property
    def protoclusters(self) -> Tuple[Protocluster, ...]:
        """ Returns the Protocluster features that the CandidateCluster was formed by """
        return tuple(self._protoclusters)

    @property
    def detection_rules(self) -> List[str]:
        """ Returns all detection rules for the child Protoluster features """
        return [cluster.detection_rule for cluster in self._protoclusters]

    @property
    def products(self) -> List[str]:
        """ Returns all unique products from contained protoclusters
            in the order they are found.
        """
        unique_products: Dict[str, None] = OrderedDict()
        for cluster in self._protoclusters:
            unique_products[cluster.product] = None
        return list(unique_products)

    @property
    def core_location(self) -> FeatureLocation:
        """ Returns a FeatureLocation covering the range of child Protocluster
            core locations
        """
        if not self._core_location:
            first_core = min(proto.core_location.start for proto in self._protoclusters)
            last_core = max(proto.core_location.end for proto in self._protoclusters)
            self._core_location = FeatureLocation(first_core, last_core)
        return self._core_location

    def get_product_string(self) -> str:
        """ Returns all unique products from contained clusters in the order
            they are found as a string, each product separated by a comma
        """
        return ",".join(self.products)

    def to_biopython(self, qualifiers: Optional[Dict[str, List[str]]] = None) -> List[SeqFeature]:
        if qualifiers is None:
            qualifiers = {}
        if self._parent_record:
            qualifiers["candidate_cluster_number"] = [str(self.get_candidate_cluster_number())]
        qualifiers["kind"] = [str(self.kind)]
        qualifiers["product"] = self.products
        qualifiers["protoclusters"] = [str(cluster.get_protocluster_number()) for cluster in self._protoclusters]
        qualifiers["detection_rules"] = self.detection_rules
        if self.smiles_structure is not None:
            qualifiers["SMILES"] = [self.smiles_structure]
        if self.polymer is not None:
            qualifiers["polymer"] = [self.polymer]
        return super().to_biopython(qualifiers)

    @classmethod
    def from_biopython(cls: Type[T], bio_feature: SeqFeature, feature: T = None,
                       leftovers: Optional[Dict] = None, record: Any = None) -> T:
        """ Does not return a proper CandidateCluster instance as extra information
            is required from the record in order to properly rebuild it
        """
        if leftovers is None:
            leftovers = Feature.make_qualifiers_copy(bio_feature)

        if not record:
            raise ValueError("record instance required for regenerating CandidateCluster from biopython")

        all_protoclusters = record.get_protoclusters()
        protocluster_numbers = [int(num) for num in leftovers.pop("protoclusters")]

        if max(protocluster_numbers) > len(all_protoclusters):
            raise ValueError("record does not contain all expected protoclusters")

        kind = CandidateClusterKind.from_string(leftovers.pop("kind")[0])
        smiles = leftovers.pop("SMILES", [None])[0]
        polymer = leftovers.pop("polymer", [None])[0]
        children = [all_protoclusters[num - 1] for num in protocluster_numbers]
        return cls(kind, children, smiles, polymer)
