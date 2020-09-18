# License: GNU Affero General Public License v3 or later
# A copy of GNU AGPL v3 should have been included in this software package in LICENSE.txt.

""" A more detailed Domain feature """

from collections import OrderedDict
from typing import Any, Dict, List, Optional, Type, TypeVar

from Bio.SeqFeature import SeqFeature

from .domain import Domain, generate_protein_location_from_qualifiers
from .feature import Feature, FeatureLocation, Location

T = TypeVar("T", bound="AntismashDomain")


class AntismashDomain(Domain):
    """ A class to represent a Domain with extra specificities and type information """
    __slots__ = ["domain_subtype", "specificity"]
    FEATURE_TYPE = "aSDomain"

    def __init__(self, location: Location, tool: str, protein_location: FeatureLocation, locus_tag: str) -> None:
        super().__init__(location, self.FEATURE_TYPE, protein_location, locus_tag,
                         tool=tool, created_by_antismash=True)
        self.domain_subtype: Optional[str] = None
        self.specificity: List[str] = []

    def to_biopython(self, qualifiers: Dict[str, List[str]] = None) -> List[SeqFeature]:
        mine: Dict[str, List[str]] = OrderedDict()
        if self.domain_subtype:
            mine["domain_subtype"] = [self.domain_subtype]
        if self.specificity:
            mine["specificity"] = self.specificity
        if qualifiers:
            mine.update(qualifiers)
        return super().to_biopython(mine)

    @classmethod
    def from_biopython(cls: Type[T], bio_feature: SeqFeature, feature: T = None,
                       leftovers: Dict[str, List[str]] = None, record: Any = None) -> T:
        if leftovers is None:
            leftovers = Feature.make_qualifiers_copy(bio_feature)
        # grab mandatory qualifiers and create the class
        tool = leftovers.pop("aSTool")[0]
        protein_location = generate_protein_location_from_qualifiers(leftovers, record)
        # locus tag is special, antismash versions <= 5.0 didn't require it, but > 5.0 do
        locus_tag = leftovers.pop("locus_tag", ["(unknown)"])[0]
        feature = cls(bio_feature.location, tool, protein_location, locus_tag)

        # grab optional qualifiers
        feature.domain_subtype = leftovers.pop("domain_subtype", [""])[0] or None
        feature.specificity = leftovers.pop("specificity", [])

        # grab parent optional qualifiers
        super().from_biopython(bio_feature, feature=feature, leftovers=leftovers, record=record)

        return feature
