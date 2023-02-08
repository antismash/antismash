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
    __slots__ = ["subtypes"]
    FEATURE_TYPE = "aSDomain"

    def __init__(self, location: Location, protein_location: FeatureLocation, locus_tag: str) -> None:
        super().__init__(location, protein_location=protein_location, locus_tag=locus_tag,
                         tool=TOOL)
        self.subtypes: List[str] = []
        self.specificity: List[str] = []

    @property
    def domain_subtype(self) -> Optional[str]:
        """ The primary subtype of the domain, if any.
            Primarily exists for backwards compatability.
        """
        if self.subtypes:
            return self.subtypes[0]
        return None

    def to_biopython(self, qualifiers: Dict[str, List[str]] = None) -> List[SeqFeature]:
        mine: Dict[str, List[str]] = OrderedDict()
        if self.subtypes:
            mine["domain_subtypes"] = self.subtypes
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
        feature.subtypes = leftovers.pop("domain_subtypes", [])
        # just in case, try for older annotation style subtype
        if not feature.subtypes and "domain_subtype" in leftovers:
            feature.subtypes.append(leftovers.pop("domain_subtype")[0])
        feature.specificity = leftovers.pop("specificity", [])

        # grab parent optional qualifiers
        super().from_biopython(bio_feature, feature=feature, leftovers=leftovers, record=record)

        return feature


register_asdomain_variant(TOOL, ModularDomain)
