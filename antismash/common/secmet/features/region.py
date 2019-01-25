# License: GNU Affero General Public License v3 or later
# A copy of GNU AGPL v3 should have been included in this software package in LICENSE.txt.

""" A class for region features """

from collections import OrderedDict
import os
from typing import Any, Dict, List, Optional, Tuple
from typing import Set  # used in comment hints, pylint: disable=unused-import
import warnings

from Bio.SeqFeature import SeqFeature
from Bio.SeqRecord import SeqRecord
from helperlibs.bio import seqio

from .cdscollection import CDSCollection, CDSFeature
from .cluster import Cluster
from .feature import Feature, FeatureLocation
from .subregion import SubRegion
from .supercluster import SuperCluster
from ..locations import combine_locations


class TemporaryRegion:
    """ A construction for the delayed conversion of a Region from biopython,
        as it requires that other feature types (SubRegion and SuperCluster)
        have already been rebuilt.

        Converts to a real Region feature with the convert_to_real_feature method.
    """
    def __init__(self, location: FeatureLocation, supercluster_numbers: List[int],
                 subregion_numbers: List[int],
                 rules: List[str] = None, products: List[str] = None,
                 probabilities: List[float] = None) -> None:
        self.type = "region"
        self.location = location
        self.superclusters = supercluster_numbers
        self.subregions = subregion_numbers
        self.detection_rules = rules
        self.products = products
        self.probabilities = probabilities

    # record should be of type Record, but that would cause a circular dependency
    def convert_to_real_feature(self, record: Any) -> "Region":
        """ Constructs a Region from this TemporaryRegion, requires the parent
            Record instance containing all the expected children of the Region
        """
        all_supers = record.get_superclusters()
        supers = [all_supers[num - 1] for num in self.superclusters]
        all_subs = record.get_subregions()
        subs = [all_subs[num - 1] for num in self.subregions]
        return Region(supers, subs)


class Region(CDSCollection):
    """ A feature that represents a region of interest made up of overlapping
        SuperCluster features and/or SubRegion features.

        At least one SuperCluster or SubRegion is required to make up a Region.

        Region features cannot overlap.
    """
    __slots__ = ["_subregions", "_superclusters", "clusterblast",
                 "knownclusterblast", "subclusterblast"]

    def __init__(self, superclusters: List[SuperCluster] = None,
                 subregions: List[SubRegion] = None) -> None:
        if not superclusters and not subregions:
            raise ValueError("A Region requires at least one child SubRegion or SuperCluster")

        children = []  # type: List[CDSCollection]

        if subregions is None:
            subregions = []
        if superclusters is None:
            superclusters = []

        for region in subregions:
            assert isinstance(region, SubRegion), type(region)
            children.append(region)
        for cluster in superclusters:
            assert isinstance(cluster, SuperCluster), type(cluster)
            children.append(cluster)

        location = combine_locations(child.location for child in children)

        super().__init__(location, feature_type="region", child_collections=children)
        self._subregions = subregions
        self._superclusters = superclusters

        self.clusterblast = None  # type: Optional[List[str]]
        self.knownclusterblast = None  # type: Any
        self.subclusterblast = None  # type: Optional[List[str]]

    @property
    def subregions(self) -> Tuple[SubRegion, ...]:
        """ Returns a list of SubRegion features used to create this region
        """
        return tuple(self._subregions)

    @property
    def superclusters(self) -> Tuple[SuperCluster, ...]:
        """ Returns a list of SuperCluster features used to create this region
        """
        return tuple(self._superclusters)

    @property
    def products(self) -> List[str]:
        """ Returns a list of unique products collected from all contained
            SuperClusters
        """
        products = OrderedDict()  # type: Dict[str, None]
        for cluster in self._superclusters:
            for product in cluster.products:
                products[product] = None
        return list(products) or ["unknown"]

    def get_product_string(self) -> str:
        """ Returns a string of all unique products collected from all
            contained SuperClusters
        """
        return "-".join(sorted(self.products))

    @property
    def detection_rules(self) -> List[str]:
        """ Returns a list of unique detection rules collected from all
            contained SuperClusters
        """
        rules = OrderedDict()  # type: Dict[str, str]
        for cluster in self._superclusters:
            for product, rule in zip(cluster.products, cluster.detection_rules):
                rules[product] = rule
        return list(rules.values())

    @property
    def probabilities(self) -> List[float]:
        """ Returns a list of probabilities collected from all contained
            SubRegions that have a probability
        """
        return [sub.probability for sub in self._subregions if sub.probability is not None]

    def add_cds(self, cds: CDSFeature) -> None:
        """ Adds a CDS to the Region and all relevant child collections. Links
            the CDS back to this region
        """
        super().add_cds(cds)
        cds.region = self

    def get_region_number(self) -> int:
        """ Returns the region's numeric ID, only guaranteed to be consistent
            when the same clusters and subregions are defined in the parent record
        """
        if not self._parent_record:
            raise ValueError("Region not in a record")
        return self._parent_record.get_region_number(self)

    def get_unique_clusters(self) -> List[Cluster]:
        """ Returns all Clusters contained by SuperClusters in this region,
            without duplicating them if multiple SuperClusters contain the same
            Cluster
        """
        clusters = set()  # type: Set[Cluster]
        for supercluster in self._superclusters:
            clusters.update(supercluster.clusters)
        return sorted(clusters)

    def write_to_genbank(self, filename: str = None, directory: str = None, record: SeqRecord = None) -> None:
        """ Writes a genbank file containing only the information contained
            within the Region.
        """
        if not filename:
            filename = "%s.region%03d.gbk" % (self.parent_record.id, self.get_region_number())
        if directory:
            filename = os.path.join(directory, filename)

        if record is None:
            record = self.parent_record.to_biopython()
        assert isinstance(record, SeqRecord)
        with warnings.catch_warnings():
            warnings.simplefilter("ignore")
            cluster_record = record[self.location.start:self.location.end]

        cluster_record.annotations["date"] = record.annotations.get("date", '')
        cluster_record.annotations["source"] = record.annotations.get("source", '')
        cluster_record.annotations["organism"] = record.annotations.get("organism", '')
        cluster_record.annotations["taxonomy"] = record.annotations.get("taxonomy", [])
        cluster_record.annotations["data_file_division"] = record.annotations.get("data_file_division", 'UNK')
        cluster_record.annotations["comment"] = record.annotations.get("comment", '')

        # update the antiSMASH annotation to include some cluster details
        comment_end_marker = "##antiSMASH-Data-END"
        cluster_comment = ("NOTE: This is a single cluster extracted from a larger record!\n"
                           "Orig. start  :: {start}\n"
                           "Orig. end    :: {end}\n"
                           "{end_marker}").format(start=self.location.start,
                                                  end=self.location.end,
                                                  end_marker=comment_end_marker)
        original = cluster_record.annotations["comment"]
        cluster_record.annotations["comment"] = original.replace(comment_end_marker, cluster_comment)

        # our cut-out clusters are always linear
        cluster_record.annotations["topology"] = "linear"

        # renumber clusters, superclusters and regions to reflect changes
        first_supercluster = min(sc.get_supercluster_number() for sc in self.superclusters)
        first_cluster = min(cluster.get_cluster_number() for cluster in self.get_unique_clusters())
        first_subregion = min(sub.get_subregion_number() for sub in self.subregions) if self.subregions else 0
        for feature in cluster_record.features:
            if feature.type == "region":
                supers = feature.qualifiers.get("supercluster_numbers")
                if not supers:
                    continue
                feature.qualifiers["supercluster_numbers"] = [str(int(num) - first_supercluster) for num in supers]
            elif feature.type == "supercluster":
                new = str(int(feature.qualifiers["supercluster_number"][0]) - first_supercluster)
                feature.qualifiers["supercluster_number"] = [new]
                new_clusters = [str(int(num) - first_cluster) for num in feature.qualifiers["child_cluster"]]
                feature.qualifiers["child_cluster"] = new_clusters
            elif feature.type in ["cluster", "cluster_core"]:
                new = str(int(feature.qualifiers["cluster_number"][0]) - first_cluster)
                feature.qualifiers["cluster_number"] = [new]
            elif feature.type == "subregion":
                new = str(int(feature.qualifiers["subregion_number"][0]) - first_subregion)
                feature.qualifiers["subregion_number"] = [new]

        seqio.write([cluster_record], filename, 'genbank')

    def to_biopython(self, qualifiers: Optional[Dict[str, List[str]]] = None) -> List[SeqFeature]:
        if not qualifiers:
            qualifiers = {}
        qualifiers["product"] = self.products
        qualifiers["rules"] = self.detection_rules
        qualifiers["probabilities"] = ["%.4f" % prob for prob in self.probabilities]
        qualifiers["subregion_numbers"] = [str(sub.get_subregion_number()) for sub in self._subregions]
        qualifiers["supercluster_numbers"] = [str(sup.get_supercluster_number()) for sup in self._superclusters]

        return super().to_biopython(qualifiers)

    @staticmethod
    def from_biopython(bio_feature: SeqFeature, feature: TemporaryRegion = None,  # type: ignore
                       leftovers: Dict[str, List[str]] = None) -> TemporaryRegion:
        if leftovers is None:
            leftovers = Feature.make_qualifiers_copy(bio_feature)
        return TemporaryRegion(bio_feature.location,
                               [int(num) for num in leftovers.pop("supercluster_numbers", [])],
                               [int(num) for num in leftovers.pop("subregion_numbers", [])],
                               leftovers["rules"],
                               leftovers["product"],
                               [float(prob) for prob in leftovers.pop("probabilities", [])],
                               )
