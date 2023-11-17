# License: GNU Affero General Public License v3 or later
# A copy of GNU AGPL v3 should have been included in this software package in LICENSE.txt.

""" Convertors for areas that include layout information which packs the areas
    into as small an area as possible.
"""

import dataclasses
import itertools
from typing import Any, Iterable, Optional, Self

from antismash.common.json import JSONCompatible
from antismash.common.secmet.features import (
    CandidateCluster,
    Feature,
    Region,
)
from antismash.common.secmet.features.protocluster import (
    Protocluster,
    SideloadedProtocluster,
)
from antismash.common.secmet.features.subregion import (
    SubRegion,
    SideloadedSubRegion,
)
from antismash.common.secmet.features.cdscollection import CDSCollection


@dataclasses.dataclass(kw_only=True, slots=True)
class Area:  # pylint: disable=too-many-instance-attributes
    """ A data structure for generic areas that will be displayed in the region view """
    start: int
    end: int
    kind: str
    height: int = 0
    neighbouring_start: int = -1
    neighbouring_end: int = -1
    product: str = ""
    prefix: str = ""
    category: str = ""
    tool: str = ""
    group: int = 0

    def __post_init__(self) -> None:
        if self.neighbouring_start == -1:
            self.neighbouring_start = self.start
        if self.neighbouring_end == -1:
            self.neighbouring_end = self.end

    @classmethod
    def from_feature(cls, feature: Feature, **kwargs: Any) -> Self:
        """ Builds an instance with the relevant information from a record feature """
        result = cls(start=int(feature.start), end=int(feature.end), kind=str(feature.FEATURE_TYPE), **kwargs)
        if isinstance(feature, Protocluster):
            result.start = int(feature.core_start)
            result.end = int(feature.core_end)
            assert result.neighbouring_start == feature.start
            assert result.neighbouring_end == feature.end
            result.product = feature.product
            result.category = feature.product_category
            result.tool = feature.tool
            if isinstance(feature, SideloadedProtocluster):
                result.prefix = f"{feature.tool}:"
        elif isinstance(feature, CandidateCluster):
            result.product = f"CC {feature.get_candidate_cluster_number()}: {feature.kind}"
            result.kind = "candidatecluster"
        elif isinstance(feature, SubRegion):
            result.product = feature.label
            if isinstance(feature, SideloadedSubRegion):
                result.prefix = f"{feature.tool}{':' if feature.label else ''}"
            else:
                result.tool = feature.tool
        # in case any of the above overwrote the extra args, set them once more
        # this is only vaguely safe because this is a dataclass
        for name, value in kwargs.items():
            setattr(result, name, value)
        return result

    def clone(self) -> "Area":
        """ Clones the instance, setting the 'group' attribute of both original and clone
            to the ID of the original.
        """
        if not self.group:
            self.group = id(self)
        return Area(**dataclasses.asdict(self))  # includes the newly modified group attribute

    def crosses_origin(self) -> bool:
        """ Returns True if the area crosses the origin """
        return self.neighbouring_start > self.neighbouring_end

    def offset(self, distance: int) -> None:
        """ Shifts all of the area's coordinates by the given distance """
        self.start += distance
        self.end += distance
        self.neighbouring_start += distance
        self.neighbouring_end += distance

    def to_minimal_json(self) -> dict[str, JSONCompatible]:
        """ Constructs a JSON compatible dictionary containing only values that have been
            explicitly set.
        """
        base = {key: val for key, val in dataclasses.asdict(self).items() if val != ""}
        if self.neighbouring_start == self.start:
            base.pop("neighbouring_start")
        if self.neighbouring_end == self.end:
            base.pop("neighbouring_end")
        if self.group == 0:
            base.pop("group")
        return base


@dataclasses.dataclass(kw_only=True)
class Row:
    """ A 1-dimensional bin for packing areas with coordinates.
        Only enough information is kept to track a single empty space.
    """
    start: int = 0
    end: int = -1
    initial_contents: dataclasses.InitVar[Optional[list[Feature]]] = None
    _contents: list[Feature] = dataclasses.field(default_factory=list, init=False)

    def __post_init__(self, initial_contents: Optional[list[Feature]] = None) -> None:
        if initial_contents:
            for area in initial_contents:
                self.add(area)

    @property
    def contents(self) -> tuple[Feature, ...]:
        """ Returns the features contained in the row """
        return tuple(self._contents)

    def can_fit(self, area: Feature) -> bool:
        """ Returns True if the given area can fit into the row """
        if not self._contents:
            return self.start <= area.start
        if area.crosses_origin():
            if self.end != -1:
                return False  # origin is already blocked
            return not area.overlaps_with(self._contents[0])
        return area.start > self.start and (self.end < 0 or area.end < self.end)

    def add(self, area: Feature) -> None:
        """ Adds the given area to the row after existing elements and updates the space available.
        """
        if not self.can_fit(area):
            raise ValueError(f"cannot fit area into row: {area}")
        self._contents.append(area)
        # remembering to strip ambiguous coordinates
        self.start = max(self.start, int(area.end) + 1)
        if area.crosses_origin():
            # the feature's start isn't the highest coordinate, so alter that
            self.start = int(area.location.end)
            self.end = int(area.start) - 1


def pack(areas: Iterable[CDSCollection], length: int = -1) -> list[Row]:
    """ Packs the given areas into as few rows as possible while being efficient.
        No single row will contains two areas that overlap.

        Arguments:
            areas: the areas pack
            length: an upper limit to bin coordinate size, not area capacity

        Returns:
            a list of rows

    """
    # this isn't quite a one-dimensional bin packing problem, as the items have
    # coordinates, but it uses something similar to the first-fit algorithm
    if not areas:
        return []
    rows = [Row(end=length)]
    for area in areas:
        for row in rows:
            if row.can_fit(area):
                row.add(area)
                break
        else:
            rows.append(Row())
            rows[-1].add(area)
    return rows


def adjust_cross_origin_area(area: Area, feature: Feature,
                             region_crosses_origin: bool, length: int,
                             ) -> Optional[Area]:
    """ Adjusts a cross-origin area to ensure that the coordinates are within the region to be
        drawn. As the parent region may not cross the origin itself (i.e. it spans the full record),
        the area may need to be split into two. If a split does occur, a new area will be created
        and returned, with both original and new will have the same group identifier for use in
        linking the two during visualisation.

        There are cases in which the existing coordinates need no adjustment, e.g. the area
        is contained fully within the pre-origin section of a region.

        Arguments:
            area: the area to be modified
            feature: the feature that the area was built from
            region_crosses_origin: whether the parent region crossed the origin
            length: the length of the record, for the purposes of extending coordinates

        Returns:
            a new area, if a split was necessary, otherwise None
    """
    if not (feature.crosses_origin() and area.crosses_origin()):
        raise ValueError("both area and feature must cross the origin")

    # the region a candidate or protocluster belongs to may not cross the origin itself, as it could
    # also cover the entire record (e.g. in smaller plasmids)
    # in those cases, the area's JSON must be split into two chunks

    extra: Optional[Area] = None

    def default() -> Optional[Area]:
        extra: Optional[Area] = None
        if region_crosses_origin:
            # extend everything but left neighbourhood
            area.start += length
            area.end += length
            area.neighbouring_end += length
        else:
            # split left neighbourhood
            extra = area.clone()
            area.start = length
            area.end = length
            area.product = ""
            area.neighbouring_end = length

            extra.neighbouring_start = 0
        return extra

    if not hasattr(feature, "core_start"):
        if region_crosses_origin:
            area.end += length
            area.neighbouring_end = area.end
            return None

        extra = area.clone()
        extra.start = extra.neighbouring_start = 0
        extra.end = extra.neighbouring_end = int(feature.end)

        area.end = length
        area.neighbouring_end = area.end

        return extra

    assert hasattr(feature, "core_start") and hasattr(feature, "core_end")

    # handle the case where the core crosses the origin
    if feature.core_start > feature.core_end:
        if region_crosses_origin:
            area.end += length
            area.neighbouring_end += length
        else:
            # split core
            extra = area.clone()
            area.end = length
            area.neighbouring_end = length
            extra.start = 0
            extra.neighbouring_start = 0
    # or the case where only the neighbourhood to the right of the core crosses the origin
    elif length - feature.core_start < feature.core_end:
        if region_crosses_origin:
            area.neighbouring_end += length
        else:
            # split right neighbourhood
            extra = area.clone()
            area.neighbouring_end = length
            extra.start = 0
            extra.end = 0
            extra.neighbouring_start = 0
            extra.product = ""
    # otherwise only the neighbourhood to the left of the core crosses the origin
    else:
        extra = default()

    return extra


def build_area_rows(region: Region, record_length: int, circular: bool = False) -> list[JSONCompatible]:
    """ Converts the areas within a region into JSON compatible forms to allow javascript
        to (relatively) easily draw the given region.

        Arguments:
            region: the region of interest
            record_length: the record length
            circular: whether the record is circular

        Returns:
            a list of JSON compatible data, one for each area
    """
    extend_over_origin = circular and (region.crosses_origin() or region.start == 0 and region.end == record_length)

    sub_rows = pack(region.subregions)
    candidates_to_include = []
    for candidate in region.candidate_clusters:
        if region.subregions or candidate.kind != candidate.kinds.SINGLE:
            candidates_to_include.append(candidate)
    cand_rows = pack(candidates_to_include)
    proto_rows = pack(region.get_unique_protoclusters())

    converted: list[Area] = []

    def add_area_from_feature(feature: Feature) -> None:
        new = Area.from_feature(feature, height=height)
        # areas crossing the origin have to be adjusted, either by splitting or
        # changing coordinates
        if extend_over_origin and feature.crosses_origin():
            assert new.crosses_origin()
            extra = adjust_cross_origin_area(new, feature, region.crosses_origin(), record_length)
            if extra:
                converted.append(new)
                new = extra
        # with a cross-origin region and an area in the post-origin portion, the
        # coordinates must also be adjusted to place the cluster correctly
        elif extend_over_origin and feature.is_contained_by(region.location.parts[-1]):
            if region.crosses_origin():
                new.offset(record_length)
        converted.append(new)

    height = 0

    for row in itertools.chain(cand_rows, sub_rows):
        for feature in row.contents:
            add_area_from_feature(feature)
        height += 2
    # in the case of neither candidates nor subregions present, bump the height by one
    # to account for the first label
    if not converted:
        height += 1
    for row in proto_rows:
        for proto in row.contents:
            add_area_from_feature(proto)
        height += 2

    return [area.to_minimal_json() for area in converted]
