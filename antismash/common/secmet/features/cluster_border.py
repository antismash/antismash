# License: GNU Affero General Public License v3 or later
# A copy of GNU AGPL v3 should have been included in this software package in LICENSE.txt.

""" A feature to represent a cluster border """

from collections import OrderedDict
from typing import Dict, List, Optional

from Bio.SeqFeature import SeqFeature

from .feature import Feature, FeatureLocation


class ClusterBorder(Feature):
    """ A feature representing a cluster border """
    __slots__ = ["tool", "probability", "cutoff", "extent", "product", "rule",
                 "contig_edge", "high_priority_product"]

    def __init__(self, location: FeatureLocation, tool: str, probability: float = None,
                 cutoff: int = 0, extent: int = 0,
                 product: Optional[str] = None, rule: Optional[str] = None,
                 contig_edge: bool = False, high_priority_product: bool = True) -> None:
        super().__init__(location, feature_type="cluster_border",
                         created_by_antismash=True)
        # required
        self.tool = str(tool)
        # args with simple defaults
        self.high_priority_product = bool(high_priority_product)
        self.contig_edge = bool(contig_edge)
        self.cutoff = int(cutoff)
        self.extent = int(extent)

        # more complicated args
        if product is not None:
            assert isinstance(product, str), type(product)
        self.product = product

        # specific to cluster finder
        self.probability = None
        if probability is not None:
            self.probability = float(probability)

        # specific to rule-based
        if rule is not None:
            assert isinstance(rule, str), type(rule)
        self.rule = rule

    def to_biopython(self, qualifiers: Dict[str, List[str]] = None) -> List[SeqFeature]:
        mine = OrderedDict()  # type: Dict[str, List[str]]
        mine["aStool"] = [self.tool]
        mine["contig_edge"] = [str(self.contig_edge)]
        if self.probability is not None:
            mine["probability"] = [str(self.probability)]
        if self.product:
            mine["product"] = [self.product]
        if self.cutoff:
            mine["cutoff"] = [str(self.cutoff)]
        if self.extent:
            mine["extent"] = [str(self.extent)]
        if self.rule:
            mine["rule"] = [self.rule]
        if qualifiers:
            mine.update(qualifiers)
        return super().to_biopython(mine)

    @staticmethod
    def from_biopython(bio_feature: SeqFeature, feature: "ClusterBorder" = None,    # type: ignore
                       leftovers: Dict[str, List[str]] = None) -> "ClusterBorder":
        if leftovers is None:
            leftovers = Feature.make_qualifiers_copy(bio_feature)

        # grab mandatory qualifiers and create the class
        tool = leftovers.pop("aStool")[0]

        # optional
        probability = float(leftovers.pop("probability")[0]) if "probability" in leftovers else None
        cutoff = int(leftovers.pop("cutoff", ["0"])[0])
        extent = int(leftovers.pop("extent", ["0"])[0])
        rule = leftovers.pop("rule", [None])[0]
        product = leftovers.pop("product", [None])[0]
        contig_edge = leftovers.pop("contig_edge", [""])[0] == "True"

        feature = ClusterBorder(bio_feature.location, tool, probability=probability,
                                cutoff=cutoff, extent=extent, rule=rule, product=product,
                                contig_edge=contig_edge)

        # grab parent optional qualifiers
        super(ClusterBorder, feature).from_biopython(bio_feature, feature=feature, leftovers=leftovers)

        return feature

    def __str__(self) -> str:
        return repr(self)

    def __repr__(self) -> str:
        return "ClusterBorder(%s, %s)" % (self.product, self.location)
