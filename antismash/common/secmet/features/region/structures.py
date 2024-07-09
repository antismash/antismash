# License: GNU Affero General Public License v3 or later
# A copy of GNU AGPL v3 should have been included in this software package in LICENSE.txt.

""" Datastructures for regions within records """

from collections import OrderedDict
import os
from typing import Any, Dict, List, Optional, Set, Tuple, Type, TypeVar, Union

from Bio.SeqFeature import SeqFeature
from Bio.SeqRecord import SeqRecord

from ..abstract import AbstractRegion
from ..cdscollection import CDSCollection, CDSFeature
from ..protocluster import Protocluster, SideloadedProtocluster
from ..feature import Feature
from ..subregion import SideloadedSubRegion, SubRegion
from ..candidate_cluster import CandidateCluster
from ...locations import (
    combine_locations,
)
from .helpers import RegionData, write_to_genbank

T = TypeVar("T", bound="Region")


class Region(AbstractRegion, CDSCollection):
    """ A feature that represents a region of interest made up of overlapping
        CandidateCluster features and/or SubRegion features.

        At least one CandidateCluster or SubRegion is required to make up a Region.

        Region features cannot overlap.
    """
    __slots__ = ["_subregions", "_candidate_clusters",
                 ]
    FEATURE_TYPE = "region"

    def __init__(self, candidate_clusters: List[CandidateCluster] = None,
                 subregions: List[SubRegion] = None) -> None:
        if not candidate_clusters and not subregions:
            raise ValueError("A Region requires at least one child SubRegion or CandidateCluster")

        children: List[CDSCollection] = []

        if subregions is None:
            subregions = []
        if candidate_clusters is None:
            candidate_clusters = []

        for region in subregions:
            assert isinstance(region, SubRegion), type(region)
            children.append(region)
        for cluster in candidate_clusters:
            assert isinstance(cluster, CandidateCluster), type(cluster)
            children.append(cluster)

        location = combine_locations(child.location for child in children)

        super().__init__(location, feature_type=self.FEATURE_TYPE, child_collections=children)
        self._subregions = subregions
        self._candidate_clusters = candidate_clusters

    @property
    def subregions(self) -> Tuple[SubRegion, ...]:
        """ Returns a list of SubRegion features used to create this region
        """
        return tuple(self._subregions)

    @property
    def candidate_clusters(self) -> Tuple[CandidateCluster, ...]:
        """ Returns a list of CandidateCluster features used to create this region
        """
        return tuple(self._candidate_clusters)

    @property
    def products(self) -> List[str]:
        """ Returns a list of unique products collected from all contained
            CandidateClusters
        """
        products: Dict[str, None] = OrderedDict()
        for cluster in self._candidate_clusters:
            for product in cluster.products:
                products[product] = None
        return list(products) or ["unknown"]

    @property
    def product_categories(self) -> Set[str]:
        """ Returns a list of unique product categories collected from all contained
            CandidateClusters
        """
        categories: Set[str] = set()
        for cluster in self._candidate_clusters:
            categories.update(cluster.product_categories)
        return categories or {"unknown"}

    def get_product_string(self) -> str:
        """ Returns a string of all unique products collected from all
            contained CandidateClusters
        """
        return ",".join(sorted(self.products))

    @property
    def detection_rules(self) -> List[str]:
        """ Returns a list of unique detection rules collected from all
            contained CandidateClusters
        """
        rules: Dict[str, str] = OrderedDict()
        for cluster in self._candidate_clusters:
            for product, rule in zip(cluster.products, cluster.detection_rules):
                rules[product] = rule
        return list(rules.values())

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

    def get_unique_protoclusters(self) -> List[Protocluster]:
        """ Returns all Protoclusters contained by CandidateClusters in this region,
            without duplicating them if multiple CandidateClusters contain the same
            Protocluster

            Result is sorted by location start, then by decreasing size, then by product
        """
        clusters: Set[Protocluster] = set()
        for candidate_cluster in self._candidate_clusters:
            clusters.update(candidate_cluster.protoclusters)
        return sorted(clusters, key=lambda x: (x.location.start, -len(x.location), x.product))

    def get_sideloaded_areas(self) -> List[Union[SideloadedProtocluster, SideloadedSubRegion]]:
        """ Returns all protoclusters and subregions that were created by
            sideloaded annotations

            Result is sorted by location start, then by decreasing size
        """
        areas: List[Union[SideloadedProtocluster, SideloadedSubRegion]] = []
        for proto in self.get_unique_protoclusters():
            if isinstance(proto, SideloadedProtocluster):
                areas.append(proto)
        for sub in self._subregions:
            if isinstance(sub, SideloadedSubRegion):
                areas.append(sub)
        return sorted(areas, key=lambda x: (x.location.start, -len(x.location), x.tool))

    def write_to_genbank(self, filename: str = None, directory: str = None, record: SeqRecord = None) -> None:
        """ Writes a genbank file containing only the information contained
            within the Region.
        """
        if not filename:
            filename = f"{self.parent_record.id}.region{self.get_region_number():03d}.gbk"
        if directory:
            filename = os.path.join(directory, filename)

        if record is None:
            record = self.parent_record.to_biopython()
        assert isinstance(record, SeqRecord)

        data = RegionData(
            start=self.location.start,
            end=self.location.end,
            candidate_clusters=self.candidate_clusters,
            subregions=self.subregions,
        )
        with open(filename, "w", encoding="utf-8") as handle:
            write_to_genbank(data, record, handle)

    def to_biopython(self, qualifiers: Optional[Dict[str, List[str]]] = None) -> List[SeqFeature]:
        if not qualifiers:
            qualifiers = {}
        if self._parent_record:
            qualifiers["region_number"] = [str(self.get_region_number())]
        qualifiers["product"] = self.products
        qualifiers["rules"] = self.detection_rules
        qualifiers["subregion_numbers"] = [str(sub.get_subregion_number()) for sub in self._subregions]
        candidates = [str(cand.get_candidate_cluster_number()) for cand in self._candidate_clusters]
        qualifiers["candidate_cluster_numbers"] = candidates

        return super().to_biopython(qualifiers)

    @classmethod
    def from_biopython(cls: Type[T], bio_feature: SeqFeature, feature: T = None,
                       leftovers: Dict[str, List[str]] = None, record: Any = None) -> T:
        if leftovers is None:
            leftovers = Feature.make_qualifiers_copy(bio_feature)

        candidate_numbers = [int(num) for num in leftovers.pop("candidate_cluster_numbers", [])]
        subregion_numbers = [int(num) for num in leftovers.pop("subregion_numbers", [])]

        if not record:
            raise ValueError("record instance required for regenerating Region from biopython")

        all_candidates = record.get_candidate_clusters()
        all_subs = record.get_subregions()

        if candidate_numbers and max(candidate_numbers) > len(all_candidates):
            raise ValueError("record does not contain all expected candidate clusters")
        if subregion_numbers and max(subregion_numbers) > len(all_subs):
            raise ValueError("record does not contain all expected subregions")

        candidates = [all_candidates[num - 1] for num in candidate_numbers]
        subs = [all_subs[num - 1] for num in subregion_numbers]

        return cls(candidates, subs)
