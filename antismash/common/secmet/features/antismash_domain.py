# License: GNU Affero General Public License v3 or later
# A copy of GNU AGPL v3 should have been included in this software package in LICENSE.txt.

""" A more detailed Domain feature """

from typing import Any, Dict, List, Type, TypeVar

from Bio.SeqFeature import SeqFeature

from .domain import Domain, generate_protein_location_from_qualifiers
from .feature import Feature, FeatureLocation, Location, pop_locus_qualifier

T = TypeVar("T", bound="AntismashDomain")


class AntismashDomain(Domain):
    """ A class to represent a Domain with extra specificities and type information """
    __slots__ = ["specificity"]
    FEATURE_TYPE = "aSDomain"

    def __init__(self, location: Location, tool: str, protein_location: FeatureLocation,
                 locus_tag: str, domain: str = None) -> None:
        super().__init__(location, self.FEATURE_TYPE, protein_location, locus_tag,
                         tool=tool, created_by_antismash=True, domain=domain)

    @classmethod
    def from_biopython(cls: Type[T], bio_feature: SeqFeature, feature: T = None,
                       leftovers: Dict[str, List[str]] = None, record: Any = None) -> T:
        if leftovers is None:
            leftovers = Feature.make_qualifiers_copy(bio_feature)
        if not feature:
            # if there is an appropriate subtype and it's not yet called, call it
            try:
                tool = leftovers["aSTool"][0]
            except KeyError:
                raise ValueError(f"{cls.FEATURE_TYPE} missing expected qualifier: 'aSTool'")
            if tool in _SUBTYPE_MAPPING:
                subtype = _SUBTYPE_MAPPING[tool]
                # this to_biopython() will call this same function again
                # so this path *must* return
                variant: AntismashDomain = subtype.from_biopython(bio_feature, feature=None,
                                                                  leftovers=leftovers, record=record)
                assert isinstance(variant, AntismashDomain)
                assert variant.type == AntismashDomain.FEATURE_TYPE
                assert isinstance(variant, cls)
                return variant

            # no subtype, so process the minimum for a default domain
            # grab mandatory qualifiers and create the class
            tool = leftovers.pop("aSTool")[0]
            protein_location = generate_protein_location_from_qualifiers(leftovers, record)
            # locus tag is special, antismash versions <= 5.0 didn't require it, but > 5.0 do
            locus_tag = pop_locus_qualifier(leftovers, allow_missing=True)
            if not locus_tag:
                raise ValueError(f"{cls.FEATURE_TYPE} missing expected qualifier: 'locus_tag'")
            assert locus_tag  # even if it's just the default "(unknown)"
            feature = cls(bio_feature.location, tool, protein_location, locus_tag)

        # for any instance, populate with the superclass info
        super().from_biopython(bio_feature, feature=feature, leftovers=leftovers, record=record)
        assert feature.domain_id
        return feature


_SUBTYPE_MAPPING: Dict[str, Type[AntismashDomain]] = {}


def register_asdomain_variant(tool: str, subtype: Type[AntismashDomain]) -> None:
    """ Register a subclass of AntismashDomain for the purposes of instantiating
        and accurate regeneration. Allows modules to define their own domain types
        without requiring all modifications to be in secmet.

        Arguments:
            tool: the name of the tool/module registering the subtype
            subtype: the type to register

        Returns:
            None
    """

    if not issubclass(subtype, AntismashDomain):
        raise TypeError(f"{subtype} is not a subclass of AntismashDomain")
    if tool in _SUBTYPE_MAPPING:
        raise ValueError(f"{tool!r} is already present as a subtype ({_SUBTYPE_MAPPING[tool]})")
    _SUBTYPE_MAPPING[tool] = subtype
