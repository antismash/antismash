# License: GNU Affero General Public License v3 or later
# A copy of GNU AGPL v3 should have been included in this software package in LICENSE.txt.

""" A more detailed Domain feature """

from collections import OrderedDict
from typing import Dict, List

from Bio.SeqFeature import SeqFeature

from .domain import Domain
from .feature import Feature, FeatureLocation


class AntismashDomain(Domain):
    """ A class to represent a Domain with extra specificities and type information """
    __slots__ = ["domain_subtype", "specificity"]

    def __init__(self, location: FeatureLocation) -> None:
        super().__init__(location, feature_type="aSDomain")
        self.domain_subtype = None  # type: str
        self.specificity = []  # type: List[str]

    def to_biopython(self, qualifiers: Dict[str, List[str]] = None) -> List[SeqFeature]:
        mine = OrderedDict()  # type: Dict[str, List[str]]
        if self.domain_subtype:
            mine["domain_subtype"] = [self.domain_subtype]
        if self.specificity:
            mine["specificity"] = self.specificity
        if self.asf:
            mine["ASF"] = self.asf.to_biopython()
        if qualifiers:
            mine.update(qualifiers)
        return super().to_biopython(mine)

    @staticmethod
    def from_biopython(bio_feature: SeqFeature, feature: "AntismashDomain" = None,  # type: ignore
                       leftovers: Dict[str, List[str]] = None) -> "AntismashDomain":
        if leftovers is None:
            leftovers = Feature.make_qualifiers_copy(bio_feature)
        # grab mandatory qualifiers and create the class
        feature = AntismashDomain(bio_feature.location)

        # grab optional qualifiers
        feature.domain_subtype = leftovers.pop("domain_subtype", [None])[0]
        feature.specificity = leftovers.pop("specificity", [])

        # grab parent optional qualifiers
        super(AntismashDomain, feature).from_biopython(bio_feature, feature=feature, leftovers=leftovers)

        return feature
