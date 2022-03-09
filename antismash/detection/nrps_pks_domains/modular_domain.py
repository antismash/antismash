# License: GNU Affero General Public License v3 or later
# A copy of GNU AGPL v3 should have been included in this software package in LICENSE.txt.

""" An AntismashDomain subclass with specificity """

from collections import OrderedDict
from typing import Any, Dict, List, Optional, Type, TypeVar

from Bio.SeqFeature import SeqFeature

from antismash.common.secmet.features.antismash_domain import (
    AntismashDomain,
    Feature,
    FeatureLocation,
    Location,
    register_asdomain_variant,
    generate_protein_location_from_qualifiers,
    pop_locus_qualifier,
)

T = TypeVar("T", bound="ModularDomain")

TOOL = "nrps_pks_domains"


class ModularDomain(AntismashDomain):
    """ A class to represent a Domain with extra specificities and type information """
    __slots__ = ["specificity"]
    FEATURE_TYPE = "aSDomain"

    def __init__(self, location: Location, protein_location: FeatureLocation, locus_tag: str) -> None:
        super().__init__(location, protein_location=protein_location, locus_tag=locus_tag,
                         tool=TOOL)
        self.domain_subtype: Optional[str] = None
        self.specificity: List[str] = []

    def to_biopython(self, qualifiers: Dict[str, List[str]] = None) -> List[SeqFeature]:
        mine: Dict[str, List[str]] = OrderedDict()
        if self.domain_subtype:
            mine["domain_subtype"] = [self.domain_subtype]
        if self.specificity:
            mine["specificity"] = list(self.specificity)
        if qualifiers:
            mine.update(qualifiers)
        assert self.domain_id
        return super().to_biopython(mine)

    @classmethod
    def from_biopython(cls: Type[T], bio_feature: SeqFeature, feature: T = None,
                       leftovers: Dict[str, List[str]] = None, record: Any = None) -> T:
        if leftovers is None:
            leftovers = Feature.make_qualifiers_copy(bio_feature)
        # grab mandatory qualifiers and create the class
        tool = leftovers.pop("aSTool")[0]
        if tool != TOOL:
            raise ValueError(f"incompatible tool type for {cls}: {tool}")
        protein_location = generate_protein_location_from_qualifiers(leftovers, record)
        locus_tag = pop_locus_qualifier(leftovers)
        assert locus_tag
        feature = cls(bio_feature.location, protein_location, locus_tag)

        # grab optional qualifiers
        feature.domain_subtype = leftovers.pop("domain_subtype", [""])[0] or None
        feature.specificity = leftovers.pop("specificity", [])

        # grab parent optional qualifiers
        super().from_biopython(bio_feature, feature=feature, leftovers=leftovers, record=record)

        return feature


register_asdomain_variant(TOOL, ModularDomain)
