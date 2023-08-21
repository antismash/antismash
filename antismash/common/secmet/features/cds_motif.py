# License: GNU Affero General Public License v3 or later
# A copy of GNU AGPL v3 should have been included in this software package in LICENSE.txt.

""" A class for CDS motif features """

from collections import OrderedDict
from typing import Any, Dict, List, Optional, Type, TypeVar, cast

from Bio.SeqFeature import SeqFeature

from .domain import Domain, generate_protein_location_from_qualifiers
from .feature import Feature, FeatureLocation, Location, pop_locus_qualifier

T = TypeVar("T", bound="CDSMotif")
ExternalT = TypeVar("ExternalT", bound="ExternalCDSMotif")


class CDSMotif(Domain):
    """ A base class for features that represent a motif within a CDSFeature """
    __slots__ = ["motif"]
    FEATURE_TYPE = "CDS_motif"

    def __init__(self, location: Location, locus_tag: str, protein_location: FeatureLocation,
                 tool: str) -> None:
        super().__init__(location, self.FEATURE_TYPE, protein_location, locus_tag,
                         tool=tool, created_by_antismash=True)

    @classmethod
    def from_biopython(cls: Type[T], bio_feature: SeqFeature, feature: T = None,
                       leftovers: Optional[Dict[str, List[str]]] = None, record: Any = None) -> T:
        if leftovers is None:
            leftovers = Feature.make_qualifiers_copy(bio_feature)
        if not feature:
            tool = leftovers.pop("aSTool", [""])[0]
            if not tool:
                return cast(T, ExternalCDSMotif.from_biopython(bio_feature, None, leftovers, record))
            protein_location = generate_protein_location_from_qualifiers(leftovers, record)
            locus_tag = pop_locus_qualifier(leftovers)
            assert locus_tag
            feature = cls(bio_feature.location, locus_tag, protein_location, tool=tool)

        updated = super().from_biopython(bio_feature, feature, leftovers, record=record)
        assert updated is feature
        assert isinstance(updated, CDSMotif)
        return updated

    def to_biopython(self, qualifiers: Dict[str, List] = None) -> List[SeqFeature]:
        mine: Dict[str, List[str]] = OrderedDict()
        if qualifiers:
            mine.update(qualifiers)
        return super().to_biopython(mine)


class ExternalCDSMotif(CDSMotif):
    """ A special case of CDSMotif annotated by external tools and not containing some expected information """
    def __init__(self, location: Location, original_qualifiers: Dict[str, List]) -> None:
        super().__init__(location, CDSMotif.FEATURE_TYPE, FeatureLocation(0, 1), "external")
        self.created_by_antismash = False
        self.original_qualifiers = original_qualifiers

    @classmethod
    def from_biopython(cls: Type[ExternalT], bio_feature: SeqFeature, feature: ExternalT = None,
                       leftovers: Optional[Dict[str, List[str]]] = None, record: Any = None) -> ExternalT:
        if leftovers is None:
            leftovers = Feature.make_qualifiers_copy(bio_feature)
        if not feature:
            feature = cls(bio_feature.location, leftovers)
        return super().from_biopython(bio_feature, feature, leftovers, record=record)

    def to_biopython(self, qualifiers: Dict[str, List] = None) -> List[SeqFeature]:
        results = super().to_biopython(qualifiers)
        for qualifier in ["aSTool", "protein_start", "protein_end", "locus_tag"]:
            results[0].qualifiers.pop(qualifier)
        # replace any removed qualifier with the original
        results[0].qualifiers.update(self.original_qualifiers)
        return results
