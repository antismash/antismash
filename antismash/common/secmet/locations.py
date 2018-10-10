# License: GNU Affero General Public License v3 or later
# A copy of GNU AGPL v3 should have been included in this software package in LICENSE.txt.

""" Helper functions for location operations """

from typing import Iterable, List, Sequence, Tuple
from typing import Optional  # in comment hints, pylint: disable=unused-import

from Bio.SeqFeature import (
    AbstractPosition,
    AfterPosition,
    BeforePosition,
    CompoundLocation,
    ExactPosition,
    FeatureLocation,
    UnknownPosition,
)

# extend FeatureLocation and CompoundLocation with bio_start / bio_end until upstream implementation
FeatureLocation.bio_start = property(lambda self: self.end if self.strand == -1 else self.start)
FeatureLocation.bio_end = property(lambda self: self.start if self.strand == -1 else self.end)
CompoundLocation.bio_start = property(lambda self: self.parts[0].bio_start)
CompoundLocation.bio_end = property(lambda self: self.parts[-1].bio_end)


def convert_protein_position_to_dna(start: int, end: int, location: FeatureLocation, codon_start: int) -> Tuple[int, int]:
    """ Convert a protein position to a nucleotide sequence position for use in generating
        new FeatureLocations from existing FeatureLocations and/or CompoundLocations.

        Arguments:
            position: the position in question, must be contained by the location
            location: the location of the related feature, for handling introns/split locations
            codon_start: the offset with which the protein coding sequence starts
        Returns:
            an int representing the calculated DNA location
    """
    if not 0 <= start < end <= len(location) // 3:
        raise ValueError("Protein positions %d and %d must be contained by %s" % (start, end, location))
    dna_start = start * 3
    dna_end = end * 3

    # typically codon_start = 1, the pythonic way
    codon_start_offset = codon_start - 1

    processed = 0
    start_found = False
    end_found = False
    for part in location.parts:
        if start_found and end_found:
            break
        if not start_found and dna_start + codon_start_offset < len(part) + processed:
            start_found = True
            if part.strand == -1:
                dna_start = part.bio_start - dna_start - codon_start_offset + processed
            else:
                dna_start = part.bio_start + dna_start + codon_start_offset - processed
        if not end_found and dna_end <= len(part) + processed:
            end_found = True
            if part.strand == -1:
                dna_end = part.bio_start - dna_end - codon_start_offset + processed
            else:
                dna_end = part.bio_start + dna_end + codon_start_offset - processed
        processed += len(part)

    assert start_found
    assert end_found

    if location.strand == -1:
        return dna_end, dna_start
    else:
        return dna_start, dna_end


def build_location_from_others(locations: Sequence[FeatureLocation]) -> FeatureLocation:
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


def location_bridges_origin(location: CompoundLocation) -> bool:
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

    for i, part in enumerate(location.parts[1:]):
        if location.strand == 1:
            if part.start <= location.parts[i].end:
                return True
        else:
            if part.start >= location.parts[i].end:
                return True
    return False


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
    lower = []  # type: List[FeatureLocation]
    upper = []  # type: List[FeatureLocation]
    if location.strand == 1:
        for part in location.parts:
            if not upper or part.start > upper[-1].end:
                upper.append(part)
            else:
                lower.append(part)
    elif location.strand == -1:
        for part in location.parts:
            if not lower or part.start < lower[-1].end:
                lower.append(part)
            else:
                upper.append(part)
    else:
        raise ValueError("Cannot separate bridged location without a valid strand")

    if not (lower and upper):
        raise ValueError("Location does not bridge origin: %s" % location)

    return lower, upper


def locations_overlap(first: FeatureLocation, second: FeatureLocation) -> bool:
    """ Returns True if the two provided FeatureLocations overlap

        Arguments:
            first: the first FeatureLocation
            second: the second FeatureLocation

        Returns:
            True if the locations overlap, otherwise False
    """
    return (first.start in second or first.end - 1 in second
            or second.start in first or second.end - 1 in first)


def location_contains_other(outer: FeatureLocation, inner: FeatureLocation) -> bool:
    """ Returns True if the first of two provided FeatureLocations contains the
        second

        Arguments:
            outer: a FeatureLocation that should contain the other
            inner: a FeatureLocation to test

        Returns:
            True if outer contains inner, otherwise False
    """
    return inner.start in outer and inner.end - 1 in outer


def location_from_string(data: str) -> FeatureLocation:
    """ Converts a string, e.g. [<1:6](-), to a FeatureLocation or CompoundLocation
    """
    def parse_position(string: str) -> AbstractPosition:
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
            strand = -1  # type: Optional[int]
        elif strand_text == '+':
            strand = 1
        elif strand_text == '?':
            strand = 0
        elif '(' not in string:
            strand = None
        else:
            raise ValueError("Cannot identify strand in location: %s" % string)

        return FeatureLocation(start, end, strand=strand)

    assert isinstance(data, str), "%s, %r" % (type(data), data)

    if '{' not in data:
        return parse_single_location(data)

    # otherwise it's a compound location
    # join{[1:6](+), [10:16](+)} -> ("join", "[1:6](+), [10:16](+)")
    operator, combined_location = data[:-1].split('{', 1)

    locations = [parse_single_location(part) for part in combined_location.split(', ')]
    return CompoundLocation(locations, operator=operator)


def combine_locations(*locations: Iterable[FeatureLocation]) -> FeatureLocation:
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
