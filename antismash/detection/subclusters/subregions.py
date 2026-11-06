# License: GNU Affero General Public License v3 or later
# A copy of GNU AGPL v3 should have been included in this software package in LICENSE.txt.

""" Construction of subregions from detected subclusters

    Sub-regions built here are added to the record during the region-formation
    step of the main pipeline, letting subclusters extend existing regions or
    form new ones. How far a detected subcluster is allowed to alter region
    boundaries depends on the mode used:
     - "create": every subcluster is kept in full, so they can also create new regions
     - "extend": overlapping subclusters are kept in full, so they can extend an
                 existing region, but subclusters without any overlap are discarded
     - "clip": only the parts overlapping an area found by another detection
               module are kept, truncated to that area, so regions can never grow
"""
from enum import StrEnum, auto
from typing import Optional, Sequence

from antismash.common.secmet import Record, SubRegion
from antismash.common.secmet.features import CDSCollection
from antismash.common.secmet.locations import (
    CompoundLocation,
    FeatureLocation,
    Location,
)

# the label given to every subregion created by this module
LABEL = "subclusters"


class SubRegionMode(StrEnum):
    """ The ways in which a detected subcluster may alter region boundaries """
    CREATE = auto()
    EXTEND = auto()
    CLIP = auto()


def gather_foreign_areas(record: Record) -> list[CDSCollection]:
    """ Returns the protoclusters and subregions of a record, all of which
        belong to other detection modules

        This module's own protoclusters are never added to the record, so
        everything found here belongs to a module that has already run.

        Arguments:
            record: the record to gather areas from

        Returns:
            a list of areas, which may be empty
    """
    areas: list[CDSCollection] = list(record.get_protoclusters())
    areas.extend(record.get_subregions())
    return areas


def build_subregions(record: Record, areas: Sequence[CDSCollection], *, tool: str,
                     mode: SubRegionMode) -> list[SubRegion]:
    """ Builds the subregion features for a set of detected subclusters

        Arguments:
            record: the record the subclusters were detected in, which must already
                    contain the areas found by any earlier detection module
            areas: the areas detected by this module, i.e. its protoclusters
            tool: the tool name to set on each resulting subregion
            mode: how far the subclusters may alter region boundaries

        Returns:
            a list of subregions, which may be empty
    """
    if mode not in SubRegionMode:
        raise ValueError(f"Unknown subcluster region mode: {mode!r}")

    subregions = [
        SubRegion(location, tool=tool, label=LABEL)
        for location, _ in record.get_potential_regions(areas)
    ]

    # "create": all subcluster predictions are used in full for subregion formation
    if mode == SubRegionMode.CREATE:
        return subregions

    # for "extend" and "clip", we want to only consider subclusters that
    # have an overlap with some other external area
    external = [location for location, _
                in record.get_potential_regions(gather_foreign_areas(record))]

    # "extend": only subclusters overlapping a foreign area are used for subregions
    if mode == SubRegionMode.EXTEND:
        filtered_areas = [area for area in areas
                          if any(area.overlaps_with(location) for location in external)]
        return [
            SubRegion(location, tool=tool, label=LABEL)
            for location, _ in record.get_potential_regions(filtered_areas)
        ]

    # "clip": only the sections shared with an area from another module are kept,
    # so the truncated results are locations rather than protoclusters

    def clip_location(location: Location, limit: Location) -> Optional[Location]:
        # intersect each part separately, as either location may cross the origin
        parts = []
        for part in location.parts:
            for limit_part in limit.parts:
                start = max(part.start, limit_part.start)
                end = min(part.end, limit_part.end)
                if end - start > 1:
                    parts.append(FeatureLocation(start, end, location.strand))

        if not parts:
            return None
        if len(parts) > 1:
            return CompoundLocation(parts)
        return parts[0]

    clipped = []
    for subregion in subregions:
        for other_location in external:
            if not subregion.overlaps_with(other_location):
                continue
            new_location = clip_location(subregion.location, other_location)
            if new_location is None:
                continue
            clipped.append(SubRegion(new_location, tool=tool, label=LABEL))

    return clipped
