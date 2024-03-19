# License: GNU Affero General Public License v3 or later
# A copy of GNU AGPL v3 should have been included in this software package in LICENSE.txt.

""" Helper functions for location operations """

import logging
from typing import Iterable, List, Optional, Sequence, Tuple, Union

from Bio.SeqFeature import (
    Position,
    AfterPosition,
    BeforePosition,
    CompoundLocation,
    ExactPosition,
    FeatureLocation,
    SeqFeature,
    UnknownPosition,
)

from .errors import SecmetInvalidInputError

Location = Union[CompoundLocation, FeatureLocation]


def convert_protein_position_to_dna(start: int, end: int, location: Location) -> Tuple[int, int]:
    """ Convert a protein position to a nucleotide sequence position for use in generating
        new FeatureLocations from existing FeatureLocations and/or CompoundLocations.

        Arguments:
            position: the position in question, must be contained by the location
            location: the location of the related feature, for handling introns/split locations

        Returns:
            an int representing the calculated DNA location
    """
    if not 0 <= start < end <= len(location) // 3:
        raise ValueError(f"Protein positions {start} and {end} must be contained by {location}")
    if location.strand == -1:
        dna_start = location.start + len(location) - end * 3
        dna_end = location.start + len(location) - start * 3
    else:
        dna_start = location.start + start * 3
        dna_end = location.start + end * 3

    # only CompoundLocations are complicated
    if not isinstance(location, CompoundLocation):
        if not location.start <= dna_start < dna_end <= location.end:
            raise ValueError(
                f"Converted coordinates {dna_start}..{dna_end} "
                f"out of bounds for location {location}"
            )
        return dna_start, dna_end

    parts = sorted(location.parts, key=lambda x: x.start)
    gap = 0
    last_end = parts[0].start
    start_found = False
    end_found = False
    for part in parts:
        if start_found and end_found:
            break
        gap += part.start - last_end
        if not start_found and dna_start + gap in part:
            start_found = True
            dna_start = dna_start + gap
        if not end_found and dna_end + gap - 1 in part:
            end_found = True
            dna_end = dna_end + gap

        last_end = part.end

    assert start_found
    assert end_found

    if not location.start <= dna_start < dna_end <= location.end:
        raise ValueError(
            f"Converted coordinates {dna_start}..{dna_end} "
            f"out of bounds for location {location}"
        )
    return dna_start, dna_end


def build_location_from_others(locations: Sequence[Location]) -> FeatureLocation:
    """ Builds a new location from non-overlapping others.
        If location boundaries are equal, they will be merged.
        If at least one provided location is a CompoundLocation or the locations
            are not continuous, the resulting location will be a CompoundLocation.

        Arguments:
            locations: a sequence of FeatureLocations to merge

        Returns:
            a FeatureLocation if the locations are continuous, otherwise a CompoundLocation
    """
    if not locations:
        raise ValueError("at least one FeatureLocation required")
    location = locations[0]
    for loc in locations[1:]:
        if loc.start == location.end:
            new_sub = FeatureLocation(location.parts[-1].start, loc.parts[0].end, location.strand)
            if len(location.parts) > 1 or len(loc.parts) > 1:
                location = CompoundLocation(location.parts[:-1] + [new_sub] + loc.parts[1:])
            else:
                location = new_sub
        else:
            location = CompoundLocation(location.parts + loc.parts)
    return location


def get_distance_between_locations(first: Location, second: Location, wrap_point: int = None) -> int:
    """ Returns the shortest distance between the two given features, crossing
        the origin if provided.

        Overlapping features are considered to have zero distance.

        Arguments:
            first: the first location
            second: the second location
            wrap_point: the point at which locations can wrap, if given

        Returns:
            the distance between the two locations
    """
    if locations_overlap(first, second):
        return 0
    offset = 0
    if wrap_point:
        offset = wrap_point
    variants = [
        abs(first.start - second.end + offset),
        abs(first.end - second.start + offset),
        abs(second.start - first.end + offset),
        abs(second.end - first.start + offset)
    ]
    distance = min(variants)
    if wrap_point:
        distance %= wrap_point
        distance = min(distance, get_distance_between_locations(first, second))
    assert distance >= 0
    return distance


def location_bridges_origin(location: Location, allow_reversing: bool = False) -> bool:
    """ Determines if a CompoundLocation would cross the origin of a record.

        Arguments:
            location: the CompoundLocation to check

        Returns:
            False if the location does not bridge the origin or if the location
            is of indeterminate strand, otherwise True
    """
    assert isinstance(location, (FeatureLocation, CompoundLocation)), type(location)

    # if it's not compound, it can't bridge at all
    if not isinstance(location, CompoundLocation):
        return False

    # invalid strands mean direction can't be determined, may need to be an error
    if location.strand not in [1, -1]:
        return False

    def check(location: Location) -> bool:
        """ Returns True if the exon ordering is invalid for the strand """
        for i, part in enumerate(location.parts[1:]):
            prev = location.parts[i]
            if location.strand == 1:
                if prev.start > part.start:
                    return True
            else:
                if prev.start < part.start:
                    return True
        return False

    if check(location):
        # due to annotations having two alternate orders for reverse strand parts:
        # 1, 2, 3, 4 and 4, 3, 2, 1, reverse the order and try again
        if allow_reversing and location.strand == -1:
            location.parts.reverse()
            if not check(location):
                logging.warning("reversed exon order for location: %s", location)
                return False
            # swap back so it will be reported as it was
            location.parts.reverse()
        return True

    return False


def _is_valid_split(lower: List[Location], upper: List[Location], strand: int) -> bool:
    """ Returns True if the results of a split are valid:
        - mutually exclusive areas covered
        - each section must be ordered correctly for the strand
    """
    if not lower or not upper:
        return False

    # check that both sections cover a mutually exclusive area
    if locations_overlap(combine_locations(lower), combine_locations(upper)):
        return False

    # check that all components in each section are correctly ordered
    for section in [upper, lower]:
        starts = [part.start for part in section]
        if sorted(starts, reverse=(strand == -1)) != starts:
            return False
    return True


def split_origin_bridging_location(location: CompoundLocation) -> Tuple[
                                                      List[FeatureLocation], List[FeatureLocation]]:
    """ Splits a CompoundLocation into two sections.
        The first contains the low-position parts (immediately after the origin
        in a forward direction), the second handles the high-position parts.

        Arguments:
            location: the CompoundLocation to split

        Returns:
            a tuple of lists, each list containing one or more FeatureLocations
    """
    lower: List[FeatureLocation] = []
    upper: List[FeatureLocation] = []
    if location.strand == 1:
        for i, part in enumerate(location.parts):
            if not upper or part.start > upper[-1].start:
                upper.append(part)
            else:
                lower.extend(location.parts[i:])
                break
    elif location.strand == -1:
        for i, part in enumerate(location.parts):
            if not lower or part.start < lower[-1].start:
                lower.append(part)
            else:
                upper.extend(location.parts[i:])
                break
    else:
        raise ValueError("Cannot separate bridged location without a valid strand")

    if not (lower and upper):
        raise ValueError(f"Location does not bridge origin: {location}")

    if not _is_valid_split(lower, upper, location.strand):
        raise ValueError(f"cannot determine correct ordering of bridged location: {location}")

    return lower, upper


def locations_overlap(first: Location, second: Location) -> bool:
    """ Returns True if the two provided FeatureLocations overlap

        Arguments:
            first: the first FeatureLocation or CompoundLocation
            second: the second FeatureLocation or CompoundLocation

        Returns:
            True if the locations overlap, otherwise False
    """
    if isinstance(first, CompoundLocation):
        return any(locations_overlap(part, second) for part in first.parts)
    if isinstance(second, CompoundLocation):
        return any(locations_overlap(first, part) for part in second.parts)
    return (first.start in second or first.end - 1 in second
            or second.start in first or second.end - 1 in first)


def location_contains_other(outer: Location, inner: Location) -> bool:
    """ Returns True if the first of two provided FeatureLocations contains the
        second

        Arguments:
            outer: a FeatureLocation or CompoundLocation that should contain the other
            inner: a FeatureLocation or CompoundLocation to test

        Returns:
            True if outer contains inner, otherwise False
    """
    if isinstance(inner, CompoundLocation):
        return all(location_contains_other(outer, part) for part in inner.parts)
    if isinstance(outer, CompoundLocation):
        return any(location_contains_other(part, inner) for part in outer.parts)
    return outer.start <= inner.start <= inner.end <= outer.end


def location_from_string(data: str) -> Location:
    """ Converts a string, e.g. [<1:6](-), to a FeatureLocation or CompoundLocation
    """
    def parse_position(string: str) -> Position:
        """ Converts a positiong from a string into a Position subclass """
        if string[0] == '<':
            return BeforePosition(int(string[1:]))
        if string[0] == '>':
            return AfterPosition(int(string[1:]))
        if string == "UnknownPosition()":
            return UnknownPosition()
        return ExactPosition(int(string))

    def parse_single_location(string: str) -> FeatureLocation:
        """ Converts a single location from a string to a FeatureLocation """
        start = parse_position(string[1:].split(':', 1)[0])  # [<1:6](-) -> <1
        end = parse_position(string.split(':', 1)[1].split(']', 1)[0])  # [<1:6](-) -> 6

        strand_text = string[-2]  # [<1:6](-) -> -
        if strand_text == '-':
            strand: Optional[int] = -1
        elif strand_text == '+':
            strand = 1
        elif strand_text == '?':
            strand = 0
        elif '(' not in string:
            strand = None
        else:
            raise ValueError(f"Cannot identify strand in location: {string}")

        return FeatureLocation(start, end, strand=strand)

    assert isinstance(data, str), f"{type(data)}, {data!r}"

    if '{' not in data:
        return parse_single_location(data)

    # otherwise it's a compound location
    # join{[1:6](+), [10:16](+)} -> ("join", "[1:6](+), [10:16](+)")
    operator, combined_location = data[:-1].split('{', 1)

    locations = [parse_single_location(part) for part in combined_location.split(', ')]
    return CompoundLocation(locations, operator=operator)


def combine_locations(*locations: Iterable[Location]) -> Location:
    """ Combines multiple FeatureLocations into a single location using the
        minimum start and maximum end. Will not create a CompoundLocation if any
        of the inputs are CompoundLocations.

        Strand will be set to None.

        Arguments:
            locations: one or more FeatureLocation instances

        Returns:
            a new FeatureLocation that will contain all provided FeatureLocations
    """
    # ensure we have a list of featureLocations
    if len(locations) == 1:
        if isinstance(locations[0], CompoundLocation):
            locs = locations[0].parts
        # it's silly to combine a single location, but don't iterate over it
        elif isinstance(locations[0], FeatureLocation):
            locs = [locations[0]]
        else:  # some kind of iterable, hopefully containing locations
            locs = list(locations[0])
    else:
        locs = list(locations)

    # build the result
    start = min(loc.start for loc in locs)
    end = max(loc.end for loc in locs)
    return FeatureLocation(start, end, strand=None)


def location_contains_overlapping_exons(location: Location) -> bool:
    """ Checks for multiple exons with the same end location, meaning they use the
        same stop codon

        Arguments:
            location: the location to check

        Returns:
            True if the location contains exons sharing a stop codon
    """
    if isinstance(location, FeatureLocation):
        return False
    if not isinstance(location, CompoundLocation):
        raise TypeError(f"expected CompoundLocation, not {type(location)}")

    return len(set(part.end for part in location.parts)) != len(location.parts)


def ensure_valid_locations(features: List[SeqFeature], can_be_circular: bool, sequence_length: int) -> None:
    """ Checks all features for valid locations, raising a ValueError if they are not

        For a location to be considered invalid, it will be one of:
            - missing a location (biopython may strip invalid locations)
            - outside the sequence provided
            - be a CDS and contain an exon with exact positions that is less than 3 bases
            - contain exons that overlap by 3 or more bases (allows for frameshifts)
            - contain exons in an order that isn't consistent with other features
                 (barring cross-origin features in records that can be circular)

        Arguments:
            features: a list of SeqFeatures (no secmet Feature should have these issues)
            can_be_circular: whether the record containing the features can be a circular genome
            sequence_length: the length of the sequence the features belong to

        Returns:
            None
    """
    for feature in features:
        # biopython drops invalid locations, so catch that first
        if feature.location is None:
            raise ValueError("one or more features with missing or invalid locations")
        # features outside the sequence cause problems with motifs and translations
        if feature.location.end > sequence_length:
            raise ValueError("feature outside record sequence: {feature.location}")
        # features with overlapping exons cause translation problems
        if location_contains_overlapping_exons(feature.location):
            raise ValueError(f"location contains overlapping exons: {feature.location}")

    # non-circular records with compound locations need to have the right part ordering
    # for translations, only really relevant for reverse strand features
    # first find what pattern has been used for locations
    standard = 0
    non_standard = 0
    for feature in features:
        if not feature.location.strand or feature.type not in ["CDS", "gene"]:
            continue

        if location_bridges_origin(feature.location):
            non_standard += 1
        else:
            standard += 1

    if can_be_circular:
        if non_standard > 2:  # allowing for a cross origin CDS and its containing gene
            raise ValueError("inconsistent exon ordering for features")
        return

    if standard and non_standard:
        raise ValueError("inconsistent exon ordering for features in non-circular record")

    if non_standard:
        for feature in features:
            if not feature.location.strand:
                continue
            if location_bridges_origin(feature.location, allow_reversing=True):
                raise ValueError(f"cannot determine correct exon ordering for location: {feature.location}")


def _adjust_location_by_offset(location: Location, offset: int) -> Location:
    """ Adjusts the given location to account for an offset (e.g. start_codon)

        Negative values are allowed, since that allows for adjusting back to
        an original location.
    """
    if offset == 0:
        return location

    assert -2 <= offset <= 2, f"invalid offset {offset}"

    def adjust_single_location(part: FeatureLocation) -> FeatureLocation:
        """ only functions on FeatureLocation """
        assert not isinstance(part, CompoundLocation)
        start = part.start
        end = part.end
        if part.strand == -1:
            end = type(end)(end + offset)
        else:
            start = type(start)(start + offset)
        return FeatureLocation(start, end, part.strand)

    if isinstance(location, CompoundLocation):
        part = location.parts[0]
        if location.strand == -1:
            assert part.end == location.end
        else:
            assert part.start == location.start
        location = CompoundLocation([adjust_single_location(part)] + location.parts[1:])
    else:
        location = adjust_single_location(location)

    return location


def frameshift_location_by_qualifier(location: Location, raw_start: Union[str, int],
                                     undo: bool = False) -> Location:
    """ Generates a new location to represent a frameshift of an existing location.
        Forward strand locations will have their start coordinate lowered.
        Reverse strand locations will have their end coordinate raised.

        Arguments:
            location: the location to shift
            start: a 1-indexed integer or string as per the genbank "codon_start" qualifier
            undo: whether to treat the frameshift as undoing a previous frameshift

        Returns:
            a new location instance of the same type
    """
    if isinstance(raw_start, str):
        try:
            codon_start = int(raw_start[0]) - 1
        except ValueError:
            raise SecmetInvalidInputError(f"invalid codon_start qualifier: {raw_start}")
    else:
        codon_start = raw_start - 1

    if not 0 <= codon_start <= 2:
        raise SecmetInvalidInputError(f"invalid codon_start qualifier: {codon_start + 1}")

    if location.strand == -1:
        codon_start *= -1

    if undo:
        codon_start *= -1

    return _adjust_location_by_offset(location, codon_start)


def offset_location(location: Location, offset: int) -> Location:
    """ Creates a new location at the given offset to the original.
        Will not loop over the origin and offsets cannot make locations negative.

        Arguments:
            location: the location to shift
            offset: the amount to offset

        Returns:
            a new location instance
    """
    parts = location.parts
    new = []
    for part in parts:
        assert part.start + offset >= 0
        new.append(FeatureLocation(part.start + offset, part.end + offset, strand=part.strand))
    if isinstance(location, CompoundLocation):
        return CompoundLocation(new, operator=location.operator)
    return new[0]


def remove_redundant_exons(location: Location) -> Location:
    """ Generates a new location that with redudant exons removed.
        Redundant exons are those that cover a location already covered by a larger exon.

        Arguments:
            location: the location to trim

        Returns:
            a new location instance, if redundant exons are found, otherwise the existing location
    """
    if len(location.parts) == 1:
        return location

    parts_by_size = sorted(location.parts, key=lambda part: part.end - part.start, reverse=True)
    parts: List[FeatureLocation] = []
    for part in parts_by_size:
        covered = False
        for existing in parts:
            if location_contains_other(existing, part):
                covered = True
                break
        if not covered:
            parts.append(part)
    if len(parts) == 1:
        return parts[0]
    return CompoundLocation([part for part in location.parts if part in parts], operator=location.operator)
