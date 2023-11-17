# License: GNU Affero General Public License v3 or later
# A copy of GNU AGPL v3 should have been included in this software package in LICENSE.txt.

""" A specific subclass of Feature for features which themselves contain
    multiple CDS features, such as clusters and candidate clusters
"""

from collections import OrderedDict
from typing import Any, Dict, List, Optional, Sequence, Tuple, Type, TypeVar, Union

from Bio.SeqFeature import SeqFeature

from .cds_feature import CDSFeature
from .feature import Feature
from ..locations import (
    CompoundLocation,
    FeatureLocation,
    Location,
    location_bridges_origin,
    split_origin_bridging_location,
)

T = TypeVar("T", bound="CDSCollection")


class CDSCollection(Feature):
    """ A specific subclass of Feature for features which themselves contain
        multiple CDS features, such as clusters and candidate clusters.

        A location is required and may extend beyond any CDSFeatures contained
        by the CDSCollection, but cannot be smaller than any CDS added.
    """
    __slots__ = ["_parent_record", "_contig_edge", "_cdses", "_children", "_parent",
                 "_cds_cache", "_cds_cache_dirty",
                 ]

    def __init__(self, location: FeatureLocation, feature_type: str,
                 child_collections: Sequence["CDSCollection"] = None) -> None:
        super().__init__(location, feature_type, created_by_antismash=True)
        if len(location.parts) > 1 and location.strand != 1:
            raise ValueError("CDS collections spanning the origin must be in the forward strand")
        self._parent_record: Any = None  # should be Record but will cause circular dependencies
        self._contig_edge = False
        self._cdses: Dict[CDSFeature, None] = OrderedDict()
        self._cds_cache: tuple[CDSFeature, ...]
        self._cds_cache_dirty: bool = True
        self._children = child_collections
        self._parent: Optional["CDSCollection"] = None
        if self._children:
            for child in self._children:
                assert isinstance(child, CDSCollection), type(child)
                child.parent = self

    def __lt__(self, other: Union[Feature, Location]) -> bool:
        """ Collections differ from other Features in that start ties are
            resolved in the opposite order, from longest to shortest
        """

        if isinstance(other, (CompoundLocation, FeatureLocation)):
            location = other
        else:
            location = other.location

        def get_comparator(loc: Location) -> tuple[int, int]:
            start = loc.start
            if location_bridges_origin(loc):
                _, head = split_origin_bridging_location(loc)
                start = min(part.start for part in head) - max(part.end for part in head)
            return (start, -len(loc))

        return get_comparator(self.location) < get_comparator(location)

    def __contains__(self, other: Any) -> bool:
        """ Returns True if the given CDSFeature or CDSCollection is one of the
            children of this collection """
        if isinstance(other, CDSFeature):
            return other in self._cdses
        if isinstance(other, CDSCollection) and self._children:
            return any(other is child or other in child for child in self._children)
        return False

    @property
    def parent_record(self) -> Any:  # again, should be Record
        """ Returns the parent Record of the collection """
        return self._parent_record

    @parent_record.setter
    def parent_record(self, record: Any) -> None:  # again, should be Record
        """ Sets the parent record to a secmet.Record instance """
        self._parent_record = record
        self._contig_edge = self.location.start == 0 or self.location.end >= len(record.seq)

    def get_root(self) -> "CDSCollection":
        """ Returns the highest level CDSCollection that either contains this
            collection or its parent
        """
        if not self._parent:
            return self
        return self._parent.get_root()

    @property
    def parent(self) -> Optional["CDSCollection"]:
        """ The parent CDSCollection of this collection or None if this has no
            parent
        """
        return self._parent

    @parent.setter
    def parent(self, parent: Optional["CDSCollection"]) -> None:
        """ Sets the parent CDSCollection, if it is not None it must fully
            contain this collection
        """
        if parent is not None:
            assert isinstance(parent, CDSCollection), type(parent)
            assert self.is_contained_by(parent)
        self._parent = parent

    @property
    def contig_edge(self) -> bool:
        """ Returns True if the collection lies on the edge of a Record.
            Requires that parent_record has been set.
        """
        if not self._parent_record:
            raise ValueError("Cannot determine if on contig edge without parent record")
        if self.crosses_origin():
            return False
        if self._contig_edge:
            return self._contig_edge
        if self._children:
            return any(child.contig_edge for child in self._children)
        return False

    def crosses_origin(self) -> bool:  # overriding base class, since these are much more controlled
        return len(self.location.parts) > 1

    def add_cds(self, cds: CDSFeature) -> None:
        """ Add a CDS to the collection covered by this feature, also adds to
            any child collections which also contain the CDS feature
        """
        if not isinstance(cds, CDSFeature):
            raise TypeError("CDS to add not a CDSFeature")
        if not cds.is_contained_by(self):
            raise ValueError("CDS added is not contained by collection")
        self._cdses[cds] = None
        self._cds_cache_dirty = True
        if not self._children:
            return
        for child in self._children:
            if cds.is_contained_by(child):
                child.add_cds(cds)

    @property
    def cds_children(self) -> Tuple[CDSFeature, ...]:
        """ Returns the CDSFeatures that have been added to this collection,
            in the order they were added """
        if self._cds_cache_dirty:
            self._cds_cache = tuple(self._cdses)
            self._cds_cache_dirty = False
        return self._cds_cache

    @classmethod
    def from_biopython(cls: Type[T], bio_feature: SeqFeature, feature: T = None,
                       leftovers: Dict[str, List[str]] = None, record: Any = None) -> T:
        assert issubclass(cls, CDSCollection)
        if leftovers is None:
            leftovers = Feature.make_qualifiers_copy(bio_feature)

        contig_edge = leftovers.pop("contig_edge", [""])[0] == "True"
        if not feature:
            feature = cls(bio_feature.location, bio_feature.type)
            feature._contig_edge = contig_edge

        # grab parent optional qualifiers
        super().from_biopython(bio_feature, feature=feature, leftovers=leftovers, record=record)
        return feature

    def to_biopython(self, qualifiers: Optional[Dict[str, List[str]]] = None) -> List[SeqFeature]:
        if not qualifiers:
            qualifiers = {}
        if self.parent_record:
            qualifiers["contig_edge"] = [str(self.contig_edge)]
        return super().to_biopython(qualifiers)
