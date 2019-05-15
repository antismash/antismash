# License: GNU Affero General Public License v3 or later
# A copy of GNU AGPL v3 should have been included in this software package in LICENSE.txt.

""" Classes and helpers for candidate cluster features

CandidateClusters contain one or more Protocluster features.

There are four kinds of candidate cluster:
    - chemical hybrid:
        contains clusters which share Cluster-defining CDSFeatures, will also
        include clusters within that shared range that do not share a CDS provided
        that they are completely contained within the candidate cluster border,
        e.g.
            ---##A###############C#---   <- Cluster 1 with definition CDSes A and C
             --##A##--                   <- Cluster 2 with definition CDS A
                      --#B#--            <- Cluster 3 with definition CDS B
                              --#C#--    <- Cluster 4 with definition CDS C
                                   -#D#- <- Cluster 5 with definition CDS D
            Since clusters 1 and 2 share a CDS that defines those clusters, a
            chemical hybrid candidate cluster exists. Clusters 1 and 4 also share a
            defining CDS, so the hybrid candidate cluster now contains clusters 1, 2 and 4.

            Cluster 3 does not share a defining CDS with either cluster 1, 2 or 4,
            but because it is interleaved into a chemical hybrid it is included
            under the assumption that it is relevant to the other clusters.
      NOTE: This may change so that there is also an 'interleaved' candidate cluster,
            where the only difference between the interleaved and the hybrid is
            that Cluster 3 would contribute it's product to the interleaved and
            not the hybrid.
    - interleaved:
        contains clusters which do not share Cluster-defining CDS features, but
        their core locations overlap,
        e.g.
            ---#A###A###A---      <- Cluster 1 with defining CDSes marked A
               ---B##B####B---    <- Cluster 2 with defining CDSes marked B
                      ---C###C--- <- Cluster 3 with defining CDSes marked C
            Since none of the clusters share any defining CDS with any other cluster,
            it is not a chemical hybrid. All three clusters would be part of an
            interleaved candidate cluster, since A overlaps with B and B overlaps with C.
    - neighbouring:
        contains clusters which transitively overlap in their neighbourhoods
        (the '-' sections in the examples above). In the chemical hybrid example,
        as all clusters overlap in some way, all 5 would be part of a neighbouring
        candidate cluster (with clusters 1-4 also being part of a hybrid candidate cluster).
        Every cluster in a 'neighbouring' cluster will also belong to one of the
        other kinds of candidate cluster.
    - single:
        the kind for all candidate clusters where only one cluster is contained,
        only exists for consistency of access. A 'single' candidate cluster will not
        exist for a cluster which is contained in either a chemical hybrid or
        an interleaved candidate cluster. In the chemical hybrid example, only cluster 5
        would be in a 'single' candidate cluster as well as in the 'neighbouring' candidate cluster

"""

from collections import OrderedDict
from enum import Enum, unique
from typing import Any, Dict, List, Optional, Tuple
from typing import Set  # comment hint, pylint: disable=unused-import

from Bio.SeqFeature import SeqFeature

from .cdscollection import CDSCollection
from .protocluster import Protocluster
from .feature import FeatureLocation, Feature
from ..locations import locations_overlap, combine_locations


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
    __slots__ = ["_protoclusters", "_kind", "smiles_structure", "polymer"]

    def __init__(self, kind: CandidateClusterKind, protoclusters: List[Protocluster],
                 smiles: str = None, polymer: str = None) -> None:
        if not protoclusters:
            raise ValueError("A CandidateCluster cannot exist without at least one Cluster")
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


def create_candidates_from_protoclusters(protoclusters: List[Protocluster]) -> List[CandidateCluster]:
    """ Constructs CandidateCluster features from a collection of Protocluster features

        CandidateClusters created may overlap if they are of different kinds.

        Chemical hybrid CandidateClusters may also contain one or more Protoclusters that
        do not share a Protocluster-defining CDS with the others, provided that the
        Protocluster(s) are fully contained within the area covered by other Protoclusters
        that do shared defining CDS features.

        Arguments:
            clusters: a list of Protoluster features

        Returns:
            a list of CandidateCluster features
    """
    if not protoclusters:
        return []

    # the implicit sort is ignored here as the cluster core location is important
    clusters = sorted(protoclusters, key=lambda x: (x.core_location.start, -len(x.core_location)))

    hybrid_area = clusters[0].core_location
    core_area = clusters[0].core_location
    neighbour_area = clusters[0].location
    # a chemical hybrid and a core overlap cannot have the same location
    # as the end point would need to use the same CDS, so a single set here is fine
    existing_locations = set()  # type: Set[Tuple[int, int]]
    hybrid_cores = set(clusters[0].definition_cdses)
    hybrid_clusters = [clusters[0]]
    core_clusters = [clusters[0]]
    neighbouring_clusters = [clusters[0]]

    # avoid creating singles for any cluster in a hybrid/overlap combo
    clusters_in_hybrids_and_overlaps = set()  # type: Set[Protocluster]

    candidate_clusters = []

    def finalise_hybrid() -> None:
        """ Construct a chemical hybrid cluster if unique location and multiple clusters """
        if len(hybrid_clusters) == 1:
            return
        candidate_cluster = CandidateCluster(CandidateClusterKind.CHEMICAL_HYBRID, hybrid_clusters)
        if (candidate_cluster.location.start, candidate_cluster.location.end) in existing_locations:
            return
        existing_locations.add((candidate_cluster.location.start, candidate_cluster.location.end))
        candidate_clusters.append(candidate_cluster)

    def finalise_nonhybrid(kind: CandidateClusterKind, clusters: List[Protocluster]) -> None:
        """ Construct a non-hybrid cluster if unique location
            and not a single cluster already in a hybrid or interleaved candidate cluster
        """
        if len(clusters) == 1:
            kind = CandidateClusterKind.SINGLE
            if clusters[0] in clusters_in_hybrids_and_overlaps:
                return
        candidate_cluster = CandidateCluster(kind, clusters)
        if (candidate_cluster.location.start, candidate_cluster.location.end) in existing_locations:
            return
        existing_locations.add((candidate_cluster.location.start, candidate_cluster.location.end))
        candidate_clusters.append(candidate_cluster)

    for cluster in clusters[1:]:
        if (locations_overlap(cluster.core_location, hybrid_area)
                and hybrid_cores.intersection(cluster.definition_cdses)):
            hybrid_area = combine_locations(cluster.core_location, hybrid_area)
            hybrid_clusters.append(cluster)
            hybrid_cores.update(cluster.definition_cdses)
            clusters_in_hybrids_and_overlaps.add(cluster)
        else:
            finalise_hybrid()
            hybrid_area = cluster.core_location
            hybrid_clusters = [cluster]
            hybrid_cores = set(cluster.definition_cdses)

        if locations_overlap(cluster.core_location, core_area):
            core_area = combine_locations(cluster.core_location, core_area)
            core_clusters.append(cluster)
            clusters_in_hybrids_and_overlaps.add(cluster)
        else:
            finalise_nonhybrid(CandidateClusterKind.INTERLEAVED, core_clusters)
            core_area = cluster.core_location
            core_clusters = [cluster]

        if locations_overlap(cluster.location, neighbour_area):
            neighbour_area = combine_locations(cluster.location, neighbour_area)
            neighbouring_clusters.append(cluster)
        else:
            finalise_nonhybrid(CandidateClusterKind.NEIGHBOURING, neighbouring_clusters)
            neighbour_area = cluster.location
            neighbouring_clusters = [cluster]

    finalise_hybrid()
    finalise_nonhybrid(CandidateClusterKind.INTERLEAVED, core_clusters)
    finalise_nonhybrid(CandidateClusterKind.NEIGHBOURING, neighbouring_clusters)

    return sorted(candidate_clusters)
