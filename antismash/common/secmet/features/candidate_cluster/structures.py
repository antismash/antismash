# License: GNU Affero General Public License v3 or later
# A copy of GNU AGPL v3 should have been included in this software package in LICENSE.txt.

""" Contains the classes required for candidate clusters """

from collections import OrderedDict
from enum import Enum, unique
from typing import Any, Dict, List, Optional, Tuple

from Bio.SeqFeature import SeqFeature

from ..cdscollection import CDSCollection
from ..protocluster import Protocluster
from ..feature import FeatureLocation, Feature
from ...locations import combine_locations


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


class TemporaryCandidateCluster:  # pylint: disable=too-many-instance-attributes
    """ A construction for the delayed conversion of a CandidateCluster from biopython,
        as it requires that other feature types (protocluster)
        have already been rebuilt.

        Converts to a real CandidateCluster feature with the convert_to_real_feature method.
    """
    def __init__(self, location: FeatureLocation, kind: CandidateClusterKind, protocluster_numbers: List[int],
                 products: List[str], detection_rules: List[str], own_number: int, contig_edge: bool,
                 smiles: str = None, polymer: str = None) -> None:  # pylint: disable=too-many-arguments
        self.type = CandidateCluster.FEATURE_TYPE
        self.location = location
        self.kind = kind
        self.protoclusters = protocluster_numbers
        self._own_number = own_number
        self.products = products
        self.detection_rules = detection_rules
        self.contig_edge = contig_edge
        self.polymer = polymer
        self.smiles_structure = smiles

    def get_candidate_cluster_number(self) -> int:
        """ Returns the candidate clusters's numeric ID, only guaranteed to be consistent for
            when the same clusters and subregions are defined in the parent record
        """
        return self._own_number

    # record type should be Record, but that ends up being a circular dependency
    def convert_to_real_feature(self, record: Any) -> "CandidateCluster":
        """ Constructs a CandidateCluster from this TemporaryCandidateCluster, requires the parent
            Record instance containing all the expected children of the CandidateCluster
        """
        if len(record.get_protoclusters()) < max(self.protoclusters):
            raise ValueError("Not all referenced clusters are present in the record")
        relevant_protoclusters = sorted([record.get_protocluster(num) for num in self.protoclusters])
        new = CandidateCluster(self.kind, relevant_protoclusters,
                               smiles=self.smiles_structure, polymer=self.polymer)
        return new


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
        unique_products = OrderedDict()  # type: Dict[str, None]
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

    @staticmethod
    def from_biopython(bio_feature: SeqFeature, feature: "CandidateCluster" = None,  # type: ignore
                       leftovers: Optional[Dict] = None) -> TemporaryCandidateCluster:
        """ Does not return a proper CandidateCluster instance as extra information
            is required from the record in order to properly rebuild it
        """
        if leftovers is None:
            leftovers = Feature.make_qualifiers_copy(bio_feature)

        kind = CandidateClusterKind.from_string(leftovers.pop("kind")[0])
        smiles = leftovers.pop("SMILES", [None])[0]
        polymer = leftovers.pop("polymer", [None])[0]
        products = leftovers.pop("product")
        own_number = leftovers.pop("candidate_cluster_number", "")
        children = [int(num) for num in leftovers.pop("protoclusters")]
        rules = leftovers.pop("detection_rules")
        edge = leftovers.pop("contig_edge", [None])[0] == "True"
        return TemporaryCandidateCluster(bio_feature.location, kind, children, products,
                                         rules, own_number, edge, smiles, polymer)
