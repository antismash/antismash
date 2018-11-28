# License: GNU Affero General Public License v3 or later
# A copy of GNU AGPL v3 should have been included in this software package in LICENSE.txt.

""" A more detailed Domain feature """

from collections import OrderedDict
from typing import Dict, List
from typing import Optional  # comment hints, pylint: disable=unused-import

from Bio.SeqFeature import SeqFeature

from .domain import Domain
from .feature import Feature, Location


class AntismashDomain(Domain):
    """ A class to represent a Domain with extra specificities and type information """
    __slots__ = ["domain_subtype", "specificity"]

    def __init__(self, location: Location, tool: str) -> None:
        super().__init__(location, feature_type="aSDomain", tool=tool, created_by_antismash=True)
        self.domain_subtype = None  # type: Optional[str]
        self.specificity = []  # type: List[str]

    def to_biopython(self, qualifiers: Dict[str, List[str]] = None) -> List[SeqFeature]:
        mine = OrderedDict()  # type: Dict[str, List[str]]
        if self.domain_subtype:
            mine["domain_subtype"] = [self.domain_subtype]
        if self.specificity:
            mine["specificity"] = self.specificity
        if qualifiers:
            mine.update(qualifiers)
        return super().to_biopython(mine)

    @staticmethod
    def from_biopython(bio_feature: SeqFeature, feature: "AntismashDomain" = None,  # type: ignore
                       leftovers: Dict[str, List[str]] = None) -> "AntismashDomain":
        if leftovers is None:
            leftovers = Feature.make_qualifiers_copy(bio_feature)
        # grab mandatory qualifiers and create the class
        tool = leftovers.pop("aSTool")[0]
        feature = AntismashDomain(bio_feature.location, tool=tool)

        # grab optional qualifiers
        feature.domain_subtype = leftovers.pop("domain_subtype", [""])[0] or None
        feature.specificity = leftovers.pop("specificity", [])

        # grab parent optional qualifiers
        super(AntismashDomain, feature).from_biopython(bio_feature, feature=feature, leftovers=leftovers)

        return feature
