# License: GNU Affero General Public License v3 or later
# A copy of GNU AGPL v3 should have been included in this software package in LICENSE.txt.

""" A specific subclass of Feature for features which themselves contain
    multiple CDS features, such as clusters and superclusters
"""

from collections import OrderedDict
from typing import Any, List, Optional, Sequence, Tuple
from typing import Dict  # used in comment hints # pylint: disable=unused-import

from Bio.SeqFeature import SeqFeature

from .feature import Feature, FeatureLocation
from .cds_feature import CDSFeature


class CDSCollection(Feature):
    """ A specific subclass of Feature for features which themselves contain
        multiple CDS features, such as clusters and superclusters.

        A location is required and may extend beyond any CDSFeatures contained
        by the CDSCollection, but cannot be smaller than any CDS added.
    """
    __slots__ = ["_parent_record", "_contig_edge", "_cdses", "_children", "_parent"]

    def __init__(self, location: FeatureLocation, feature_type: str,
                 child_collections: Sequence["CDSCollection"] = None) -> None:
        super().__init__(location, feature_type, created_by_antismash=True)
        self._parent_record = None  # type: Any  # should be Record but will cause circular dependencies
        self._contig_edge = False
        self._cdses = OrderedDict()  # type: Dict[CDSFeature, None]
        self._children = child_collections
        self._parent = None  # type: Optional["CDSCollection"]
        if self._children:
            for child in self._children:
                assert isinstance(child, CDSCollection), type(child)
                child.parent = self

    def __lt__(self, other: "Feature") -> bool:
        """ Collections differ from other Features in that start ties are
            resolved in the opposite order, from longest to shortest
        """
        if self.location.start < other.location.start:
            return True
        if self.location.start > other.location.start:
            return False
        # when starts are equal, sort by largest collection first
        return self.location.end > other.location.end

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
        return self._contig_edge

    def add_cds(self, cds: CDSFeature) -> None:
        """ Add a CDS to the collection covered by this feature, also adds to
            any child collections which also contain the CDS feature
        """
        if not isinstance(cds, CDSFeature):
            raise TypeError("CDS to add not a CDSFeature")
        if not cds.is_contained_by(self):
            raise ValueError("CDS added is not contained by collection")
        self._cdses[cds] = None
        if not self._children:
            return
        for child in self._children:
            if cds.is_contained_by(child):
                child.add_cds(cds)

    @property
    def cds_children(self) -> Tuple[CDSFeature, ...]:
        """ Returns the CDSFeatures that have been added to this collection,
            in the order they were added """
        return tuple(self._cdses)

    @staticmethod
    def from_biopython(bio_feature: SeqFeature, feature: "CDSCollection" = None,  # type: ignore
                       leftovers: Optional[Dict] = None) -> "CDSCollection":
        if leftovers is None:
            leftovers = Feature.make_qualifiers_copy(bio_feature)

        contig_edge = leftovers.pop("contig_edge", [None])[0] == "True"
        if not feature:
            feature = CDSCollection(bio_feature.location, bio_feature.type)
            feature._contig_edge = contig_edge  # pylint: disable=protected-access

        # grab parent optional qualifiers
        super(CDSCollection, feature).from_biopython(bio_feature, feature=feature, leftovers=leftovers)
        return feature

    def to_biopython(self, qualifiers: Optional[Dict[str, List[str]]] = None) -> List[SeqFeature]:
        if not qualifiers:
            qualifiers = {}
        if self.parent_record:
            qualifiers["contig_edge"] = [str(self.contig_edge)]
        return super().to_biopython(qualifiers)
