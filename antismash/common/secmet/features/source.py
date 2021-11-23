# License: GNU Affero General Public License v3 or later
# A copy of GNU AGPL v3 should have been included in this software package in LICENSE.txt.

""" A feature to represent a gene """

from typing import Any, Dict, List, Optional, Type, TypeVar

from Bio.SeqFeature import SeqFeature

from .feature import Feature, Location

T = TypeVar("T", bound="Source")


class Source(Feature):
    """ A feature for `source` features within records.
        Typically only one of these is present, but when multiple are present
        then it is meaningful.
    """
    __slots__: List[str] = []
    FEATURE_TYPE = "source"

    def __init__(self, location: Location, created_by_antismash: bool = False,
                 qualifiers: Optional[Dict[str, List[str]]] = None) -> None:
        super().__init__(location, feature_type=self.FEATURE_TYPE,
                         created_by_antismash=created_by_antismash)
        if qualifiers:
            assert isinstance(qualifiers, dict)
            self._qualifiers.update(qualifiers)

    def to_biopython(self, qualifiers: Dict[str, Any] = None) -> SeqFeature:
        """ Construct a matching SeqFeature for this Source feature """
        if not qualifiers:
            qualifiers = {}
        return super().to_biopython(qualifiers)

    @classmethod
    def from_biopython(cls: Type[T], bio_feature: SeqFeature, feature: T = None,
                       leftovers: Dict[str, List[str]] = None, record: Any = None) -> T:
        if leftovers is None:
            leftovers = Feature.make_qualifiers_copy(bio_feature)
        feature = cls(bio_feature.location)
        super().from_biopython(bio_feature, feature=feature, leftovers=leftovers, record=record)
        return feature
