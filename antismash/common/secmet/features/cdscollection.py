# License: GNU Affero General Public License v3 or later
# A copy of GNU AGPL v3 should have been included in this software package in LICENSE.txt.

""" A specific subclass of Feature for features which themselves contain
    multiple CDS features, such as clusters and candidate clusters
"""

from collections import OrderedDict
from dataclasses import dataclass, field
from enum import Enum, auto, unique
from typing import (
    Any,
    Dict,
    Iterator,
    List,
    Optional,
    Protocol,
    Sequence,
    SupportsIndex,
    Type,
    TypeVar,
    Union,
)

from Bio.SeqFeature import SeqFeature

from .cds_feature import CDSFeature
from .feature import Feature
from ..locations import (
    CompoundLocation,
    FeatureLocation,
    Location,
    Position,
    location_bridges_origin,
    split_origin_bridging_location,
)

T = TypeVar("T", bound="CDSCollection")


@unique
class CollectionSection(Enum):
    """ Defined subsections of a collection
    """
    PRE_ORIGIN = auto()
    CROSS_ORIGIN = auto()
    POST_ORIGIN = auto()


@dataclass
class _CDSCache:
    # for the purposes of deduplication while keeping order, use an ordered dict
    _features: dict[CDSFeature, None] = field(default_factory=OrderedDict)
    _cached: tuple[CDSFeature, ...] = field(default_factory=tuple)
    _dirty: bool = True

    @property
    def features(self) -> tuple[CDSFeature, ...]:
        """ The features within the cache, updating if modified since last use """
        if self._dirty:
            self._regen_cache()
            self._dirty = False
        return self._cached

    def __contains__(self, other: CDSFeature) -> bool:
        """ Returns True if the given CDSFeature is a member of this collection """
        return other in self._features

    def add_cds(self, cds: CDSFeature, _section: CollectionSection = None) -> None:
        """ Adds the given CDS to the the cache, ignoring any specified section,
            as this level of cache is not sectioned
        """
        self._features[cds] = None
        self._dirty = True

    def clear(self) -> None:
        """ Clears the cache of the features within it """
        self._features.clear()
        self._dirty = True

    def _regen_cache(self) -> None:
        """ An overridable method for regenerating the cached object itself """
        self._cached = tuple(self._features)


class _SectionedCDSTuple(tuple):
    def __new__(cls, lookup: dict[CDSFeature, None],
                pre_origin: tuple[CDSFeature, ...],
                cross_origin: tuple[CDSFeature, ...],
                post_origin: tuple[CDSFeature, ...],
                ) -> "_SectionedCDSTuple":
        return tuple.__new__(cls, (lookup, pre_origin, cross_origin, post_origin, False))

    @property
    def _lookup(self) -> dict[CDSFeature, None]:
        return super().__getitem__(0)

    @property
    def pre_origin(self) -> tuple[CDSFeature, ...]:
        """ Contains the CDS features within the pre-origin section of the area """
        return super().__getitem__(1)

    @property
    def cross_origin(self) -> tuple[CDSFeature, ...]:
        """ Contains the CDS features within the cross-origin section of the area """
        return super().__getitem__(2)

    @property
    def post_origin(self) -> tuple[CDSFeature, ...]:
        """ Contains the CDS features within the post-origin section of the area """
        return super().__getitem__(3)

    def index(self, item: Any, start: SupportsIndex = 0, stop: SupportsIndex = None) -> int:
        if stop is None:
            stop = len(self)
        # if the limits are ints, short circuits are possible
        if isinstance(start, int) and isinstance(stop, int):
            for i, cds in enumerate(self._lookup):
                if start <= i <= stop and item == cds:
                    return i
            raise IndexError("tuple index out of range")
        # if not, the slow version is required
        return tuple(self._lookup).index(item, start, stop)

    def __contains__(self, item: object) -> bool:
        return item in self._lookup

    def __len__(self) -> int:
        return len(self._lookup)

    def __iter__(self) -> Iterator[CDSFeature]:
        for cds in self._lookup:
            yield cds

    def __getitem__(self, index: Any) -> Any:
        return list(self._lookup)[index]

    def __str__(self) -> str:
        return str(self._lookup)

    def __eq__(self, other: Any) -> bool:
        return other == tuple(self._lookup)

    def __reduce__(self) -> tuple:
        # since this subclasses tuples, pickling isn't clever enough to handle the extra arguments
        # so this function is required to make the reconstruction function with the right args
        return (_SectionedCDSTuple, (self._lookup, self.pre_origin, self.cross_origin, self.post_origin))


@dataclass
class _SectionedCDSCache(_CDSCache):
    _pre_origin: _CDSCache = field(default_factory=_CDSCache)
    _cross_origin: _CDSCache = field(default_factory=_CDSCache)
    _post_origin: _CDSCache = field(default_factory=_CDSCache)

    def add_cds(self, cds: CDSFeature, section: Optional[CollectionSection] = CollectionSection.POST_ORIGIN) -> None:
        super().add_cds(cds)
        if section is CollectionSection.PRE_ORIGIN:
            self._pre_origin.add_cds(cds)
        elif section is CollectionSection.CROSS_ORIGIN:
            self._cross_origin.add_cds(cds)
        elif section is CollectionSection.POST_ORIGIN:
            self._post_origin.add_cds(cds)
        else:
            raise ValueError(f"unknown section value: {section}")

    def clear(self) -> None:
        super().clear()
        for section in (self._pre_origin, self._cross_origin, self._post_origin):
            section.clear()

    def _regen_cache(self) -> None:
        self._cached = _SectionedCDSTuple(self._features,
                                          pre_origin=self._pre_origin.features,
                                          cross_origin=self._cross_origin.features,
                                          post_origin=self._post_origin.features,
                                          )
        assert isinstance(self._cached, _SectionedCDSTuple)

    def __str__(self) -> str:
        return f"_SectionedCDSCache: {self.features}"

    def __repr__(self) -> str:
        return f"_SectionedCDSCache: {self.features}"


class CDSCollection(Feature):
    """ A specific subclass of Feature for features which themselves contain
        multiple CDS features, such as clusters and candidate clusters.

        A location is required and may extend beyond any CDSFeatures contained
        by the CDSCollection, but cannot be smaller than any CDS added.
    """
    __slots__ = ["_parent_record", "_contig_edge", "_cdses", "_children", "_parent",
                 ]

    def __init__(self, location: Location, feature_type: str,
                 child_collections: Sequence["CDSCollection"] = None) -> None:
        # it's fine to have two parts if crossing the origin
        if len(location.parts) > 1:
            # but it must only be the two halves, more indicates a problem of some kind
            assert len(location.parts) == 2, location
            # and the second half must be the one that starts at the origin itself
            if location.parts[1].start != 0:
                raise ValueError("Collections cannot have compound locations without crossing the origin")
        assert len(set(part.strand for part in location.parts)) == 1, f"mixed strands for collection: {location=}"
        super().__init__(location, feature_type, created_by_antismash=True)
        if len(location.parts) > 1 and location.strand != 1:
            raise ValueError("CDS collections spanning the origin must be in the forward strand")
        self._parent_record: Any = None  # should be Record but will cause circular dependencies
        self._contig_edge = False
        self._cdses = _SectionedCDSCache()
        self._children = child_collections or []
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
            # shortcut here if the feature is a child of this collection
            if other in self:
                return True
            location = other.location

        # if the other location is contained within this feature, a shortcut can be taken
        if self.location.contains(location) and not location.contains(self.location):
            return True

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
        if self._children:
            for child in self._children:
                if not child.parent_record:
                    child.parent_record = record

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

    def add_cds(self, cds: CDSFeature, section: CollectionSection = None) -> None:
        """ Add a CDS to the collection covered by this feature, also adds to
            any child collections which also contain the CDS feature
        """
        if not isinstance(cds, CDSFeature):
            raise TypeError("CDS to add not a CDSFeature")
        if not cds.is_contained_by(self):
            raise ValueError("CDS added is not contained by collection")

        if section is None and self.crosses_origin():
            section = CollectionSection.PRE_ORIGIN
            if cds.crosses_origin():
                section = CollectionSection.CROSS_ORIGIN
            elif cds.is_contained_by(self.location.parts[1]):
                section = CollectionSection.POST_ORIGIN

        kwargs = {}  # allow fallbacks per child type, if a section is even needed
        if section:
            kwargs["section"] = section
        self._cdses.add_cds(cds, **kwargs)

        for child in self._children:
            if cds.is_contained_by(child):
                child.add_cds(cds, **kwargs)

    @property
    def cds_children(self) -> _SectionedCDSTuple:
        """ Returns the CDSFeatures that have been added to this collection,
            in the order they were added """
        result = self._cdses.features
        assert isinstance(result, _SectionedCDSTuple)
        return result

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


# this usage of Protocol allows subclasses to define the core location arbitrarily
# rather than just as an attribute
class _CoreProtocol(Protocol):
    @property
    def core_location(self) -> Location:
        """ The location covering the collection's core genes """


class CoredCollectionMixin(_CoreProtocol):
    """ Helpers for working with core locations in those collections that define or require them
    """

    @property
    def core_end(self) -> Position:
        """ The end coordinate of the core location.

            NOTE: differs from the location.end, as that is the maximum coordinate
        """
        if self.core_location.strand != -1:
            return self.core_location.parts[-1].end
        return self.core_location.parts[0].end

    @property
    def core_start(self) -> Position:
        """ The start coordinate of the core location.

            NOTE: differs from location.start, as that is the minimum coordinate
        """
        if self.core_location.strand != -1:
            return self.core_location.parts[0].start
        return self.core_location.parts[-1].start
