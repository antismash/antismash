# License: GNU Affero General Public License v3 or later
# A copy of GNU AGPL v3 should have been included in this software package in LICENSE.txt.

""" Classes and helpers for supercluster features

Superclusters contain one or more Cluster features.

There are four kinds of supercluster:
    - chemical hybrid:
        contains clusters which share Cluster-defining CDSFeatures, will also
        include clusters within that shared range that do not share a CDS provided
        that they are completely contained within the supercluster border,
        e.g.
            ---##A###############C#---   <- Cluster 1 with definition CDSes A and C
             --##A##--                   <- Cluster 2 with definition CDS A
                      --#B#--            <- Cluster 3 with definition CDS B
                              --#C#--    <- Cluster 4 with definition CDS C
                                   -#D#- <- Cluster 5 with definition CDS D
            Since clusters 1 and 2 share a CDS that defines those clusters, a
            chemical hybrid supercluster exists. Clusters 1 and 4 also share a
            defining CDS, so the hybrid supercluster now contains clusters 1, 2 and 4.

            Cluster 3 does not share a defining CDS with either cluster 1, 2 or 4,
            but because it is interleaved into a chemical hybrid it is included
            under the assumption that it is relevant to the other clusters.
      NOTE: This may change so that there is also an 'interleaved' supercluster,
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
            interleaved supercluster, since A overlaps with B and B overlaps with C.
    - neighbouring:
        contains clusters which transitively overlap in their neighbourhoods
        (the '-' sections in the examples above). In the chemical hybrid example,
        as all clusters overlap in some way, all 5 would be part of a neighbouring
        supercluster (with clusters 1-4 also being part of a hybrid supercluster).
        Every cluster in a 'neighbouring' cluster will also belong to one of the
        other kinds of supercluster.
    - single:
        the kind for all superclusters where only one cluster is contained,
        only exists for consistency of access. A 'single' supercluster will not
        exist for a cluster which is contained in either a chemical hybrid or
        an interleaved supercluster. In the chemical hybrid example, only cluster 5
        would be in a 'single' supercluster as well as in the 'neighbouring' supercluster

"""

from collections import OrderedDict
from enum import Enum, unique
from typing import Any, Dict, List, Optional, Tuple
from typing import Set  # comment hint, pylint: disable=unused-import

from Bio.SeqFeature import SeqFeature

from .cdscollection import CDSCollection
from .cluster import Cluster
from .feature import FeatureLocation, Feature
from ..locations import locations_overlap, combine_locations


@unique
class SuperClusterKind(Enum):
    """ An Enum representing the kind of a SuperCluster.
        Allows for more flexible conversion and more robust value constraints.
    """
    SINGLE = 0
    INTERLEAVED = 1
    NEIGHBOURING = 2
    CHEMICAL_HYBRID = 3

    def __str__(self) -> str:
        return str(self.name).lower()

    @staticmethod
    def from_string(label: str) -> "SuperClusterKind":
        """ Converts a string to a SuperClusterKind instance when possible.
            Raises an error if not possible.
        """
        for value in SuperClusterKind:
            if str(value) == label:
                return value
        raise ValueError("unknown supercluster kind: %s" % label)


class TemporarySuperCluster:  # pylint: disable=too-many-instance-attributes
    """ A construction for the delayed conversion of a SuperCluster from biopython,
        as it requires that other feature types (Cluster)
        have already been rebuilt.

        Converts to a real SuperCluster feature with the convert_to_real_feature method.
    """
    def __init__(self, location: FeatureLocation, kind: SuperClusterKind, cluster_numbers: List[int],
                 products: List[str], detection_rules: List[str], own_number: int, contig_edge: bool,
                 smiles: str = None, polymer: str = None) -> None:  # pylint: disable=too-many-arguments
        self.type = "supercluster"
        self.location = location
        self.kind = kind
        self.clusters = cluster_numbers
        self._own_number = own_number
        self.products = products
        self.detection_rules = detection_rules
        self.contig_edge = contig_edge
        self.polymer = polymer
        self.smiles_structure = smiles

    def get_supercluster_number(self) -> int:
        """ Returns the superclusters's numeric ID, only guaranteed to be consistent for
            when the same clusters and subregions are defined in the parent record
        """
        return self._own_number

    # record type should be Record, but that ends up being a circular dependency
    def convert_to_real_feature(self, record: Any) -> "SuperCluster":
        """ Constructs a SuperCluster from this TemporarySuperCluster, requires the parent
            Record instance containing all the expected children of the SuperCluster
        """
        if len(record.get_clusters()) < max(self.clusters):
            raise ValueError("Not all referenced clusters are present in the record")
        relevant_clusters = sorted([record.get_cluster(num) for num in self.clusters])
        new = SuperCluster(self.kind, relevant_clusters,
                           smiles=self.smiles_structure, polymer=self.polymer)
        return new


class SuperCluster(CDSCollection):
    """ A class representing a collection of overlapping Cluster features.
        The location of a SuperCluster is defined as the minimum area that would
        contain all of the child Clusters.
    """
    kinds = SuperClusterKind
    __slots__ = ["_clusters", "_kind", "smiles_structure", "polymer"]

    def __init__(self, kind: SuperClusterKind, clusters: List[Cluster],
                 smiles: str = None, polymer: str = None) -> None:
        if not clusters:
            raise ValueError("A SuperCluster cannot exist without at least one Cluster")
        for cluster in clusters:
            assert isinstance(cluster, Cluster), type(cluster)
        if not isinstance(kind, SuperClusterKind):
            raise TypeError("argument 1 should be SuperClusterKind, had %s" % type(kind))
        location = combine_locations(cluster.location for cluster in clusters)
        super().__init__(location, feature_type="supercluster", child_collections=clusters)
        self._clusters = clusters
        self._kind = kind
        self.smiles_structure = smiles
        self.polymer = polymer

    def __repr__(self) -> str:
        return "SuperCluster(%s, %s)" % (self.location, self.kind)

    def get_supercluster_number(self) -> int:
        """ Returns the superclusters's numeric ID, only guaranteed to be consistent for
            when the same clusters and subregions are defined in the parent record
        """
        if not self._parent_record:
            raise ValueError("SuperCluster not contained in record")
        return self._parent_record.get_supercluster_number(self)

    @property
    def kind(self) -> SuperClusterKind:
        """ The kind of SuperCluster """
        return self._kind

    @property
    def clusters(self) -> Tuple[Cluster, ...]:
        """ Returns the Cluster features that the SuperCluster was formed by """
        return tuple(self._clusters)

    @property
    def detection_rules(self) -> List[str]:
        """ Returns all detection rules for the child Cluster features """
        return [cluster.detection_rule for cluster in self._clusters]

    @property
    def products(self) -> List[str]:
        """ Returns all unique products from contained clusters
            in the order they are found.
        """
        unique_products = OrderedDict()  # type: Dict[str, None]
        for cluster in self._clusters:
            unique_products[cluster.product] = None
        return list(unique_products)

    def to_biopython(self, qualifiers: Optional[Dict[str, List[str]]] = None) -> List[SeqFeature]:
        if qualifiers is None:
            qualifiers = {}
        if self._parent_record:
            qualifiers["supercluster_number"] = [str(self.get_supercluster_number())]
        qualifiers["kind"] = [str(self.kind)]
        qualifiers["product"] = self.products
        qualifiers["child_cluster"] = [str(cluster.get_cluster_number()) for cluster in self._clusters]
        qualifiers["detection_rules"] = [cluster.detection_rule for cluster in self._clusters]
        if self.smiles_structure is not None:
            qualifiers["SMILES"] = [self.smiles_structure]
        if self.polymer is not None:
            qualifiers["polymer"] = [self.polymer]
        return super().to_biopython(qualifiers)

    @staticmethod
    def from_biopython(bio_feature: SeqFeature, feature: "SuperCluster" = None,  # type: ignore
                       leftovers: Optional[Dict] = None) -> TemporarySuperCluster:
        """ Does not return a proper SuperCluster instance as extra information
            is required from the record in order to properly rebuild it
        """
        if leftovers is None:
            leftovers = Feature.make_qualifiers_copy(bio_feature)

        kind = SuperClusterKind.from_string(leftovers.pop("kind")[0])
        smiles = leftovers.pop("SMILES", [None])[0]
        polymer = leftovers.pop("polymer", [None])[0]
        products = leftovers.pop("product")
        own_number = leftovers.pop("supercluster_number", "")
        children = [int(num) for num in leftovers.pop("child_cluster")]
        rules = leftovers.pop("detection_rules")
        edge = leftovers.pop("contig_edge", [None])[0] == "True"
        return TemporarySuperCluster(bio_feature.location, kind, children, products,
                                     rules, own_number, edge, smiles, polymer)


def create_superclusters_from_clusters(clusters: List[Cluster]) -> List[SuperCluster]:
    """ Constructs SuperCluster features from a collection of Cluster features

        SuperClusters created may overlap if they are of different kinds.

        Chemical hybrid SuperClusters may also contain one or more Clusters that
        do not share a Cluster-defining CDS with the others, provided that the
        Cluster(s) are fully contained within the area covered by other Clusters
        that do shared defining CDS features.

        Arguments:
            clusters: a list of Cluster features

        Returns:
            a list of SuperCluster features
    """
    if not clusters:
        return []

    # the implicit sort is ignored here as the cluster core location is important
    clusters = sorted(clusters, key=lambda x: (x.core_location.start, -len(x.core_location)))

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
    clusters_in_hybrids_and_overlaps = set()  # type: Set[Cluster]

    superclusters = []

    def finalise_hybrid() -> None:
        """ Construct a chemical hybrid cluster if unique location and multiple clusters """
        if len(hybrid_clusters) == 1:
            return
        supercluster = SuperCluster(SuperClusterKind.CHEMICAL_HYBRID, hybrid_clusters)
        if (supercluster.location.start, supercluster.location.end) in existing_locations:
            return
        existing_locations.add((supercluster.location.start, supercluster.location.end))
        superclusters.append(supercluster)

    def finalise_nonhybrid(kind: SuperClusterKind, clusters: List[Cluster]) -> None:
        """ Construct a non-hybrid cluster if unique location
            and not a single cluster already in a hybrid or interleaved supercluster
        """
        if len(clusters) == 1:
            kind = SuperClusterKind.SINGLE
            if clusters[0] in clusters_in_hybrids_and_overlaps:
                return
        supercluster = SuperCluster(kind, clusters)
        if (supercluster.location.start, supercluster.location.end) in existing_locations:
            return
        existing_locations.add((supercluster.location.start, supercluster.location.end))
        superclusters.append(supercluster)

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
            finalise_nonhybrid(SuperClusterKind.INTERLEAVED, core_clusters)
            core_area = cluster.core_location
            core_clusters = [cluster]

        if locations_overlap(cluster.location, neighbour_area):
            neighbour_area = combine_locations(cluster.location, neighbour_area)
            neighbouring_clusters.append(cluster)
        else:
            finalise_nonhybrid(SuperClusterKind.NEIGHBOURING, neighbouring_clusters)
            neighbour_area = cluster.location
            neighbouring_clusters = [cluster]

    finalise_hybrid()
    finalise_nonhybrid(SuperClusterKind.INTERLEAVED, core_clusters)
    finalise_nonhybrid(SuperClusterKind.NEIGHBOURING, neighbouring_clusters)

    return sorted(superclusters)
