# License: GNU Affero General Public License v3 or later
# A copy of GNU AGPL v3 should have been included in this software package in LICENSE.txt.

""" Helper functions for location operations """

from typing import Iterable, List, Sequence, Tuple

from Bio.SeqRecord import (
    SeqRecord,
)

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
FeatureLocation.bio_start = property(lambda self: self.end if self.strand == -1 else self.start )
FeatureLocation.bio_end = property(lambda self: self.start if self.strand == -1 else self.end )
CompoundLocation.bio_start = property(lambda self: self.parts[0].bio_start )
CompoundLocation.bio_end = property(lambda self: self.parts[-1].bio_end )


def convert_protein_position_to_dna(start: int, end: int, location: FeatureLocation) -> Tuple[int, int]:
    """ Convert a protein position to a nucleotide sequence position for use in generating
        new FeatureLocations from existing FeatureLocations and/or CompoundLocations.

        Arguments:
            position: the position in question, must be contained by the location
            location: the location of the related feature, for handling introns/split locations

        Returns:
            an int representing the calculated DNA location
    """
    if not 0 <= start < end <= len(location) // 3:
        raise ValueError("Protein positions %d and %d must be contained by %s" % (start, end, location))
    dna_start = start * 3
    dna_end = end * 3

    processed = 0
    start_found = False
    end_found = False
    for part in location.parts:
        if start_found and end_found:
            break
        if not start_found and dna_start < len(part) + processed:
            start_found = True
            dna_start = part.bio_start - dna_start + processed if part.strand == -1 else part.bio_start + dna_start - processed
        if not end_found and dna_end <= len(part) + processed:
            end_found = True
            dna_end = part.bio_start - dna_end + processed if part.strand == -1 else part.bio_start + dna_end - processed
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


def locations_overlap(first: FeatureLocation, second: FeatureLocation) -> bool:
    """ Returns True if the two provided FeatureLocations overlap

        Arguments:
            first: the first FeatureLocation
            second: the second FeatureLocation

        Returns:
            True if the locations overlap, otherwise False
    """
    for primary in first.parts:
        for secondary in second.parts:
            # -1 to account for the non-inclusive end
            if primary.start in secondary or primary.end - 1 in secondary or secondary.start in primary or secondary.end - 1 in primary:
                return True
    # none of sublocations overlap
    return False

def location_contains_other(outer: FeatureLocation, inner: FeatureLocation) -> bool:
    """ Returns True if the first of two provided FeatureLocations contains the
        second

        Arguments:
            outer: a FeatureLocation that should contain the other
            inner: a FeatureLocation to test

        Returns:
            True if outer contains inner, otherwise False
    """
    sublocations_found = 0
    for inner_sublocation in inner.parts:
        for outer_sublocation in outer.parts:
            # -1 to account for the non-inclusive end
            if inner_sublocation.start in outer_sublocation and inner_sublocation.end - 1 in outer_sublocation:
                sublocations_found += 1
                break # break in order to avoid scoring the query sublocation twice in weirdly overlapping subjects
    return len(inner.parts) == sublocations_found


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
            strand = -1
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


def extend_location_by(location: FeatureLocation, by: int, record: SeqRecord) -> FeatureLocation:
    """ Extends a single FeatureLocation by N bp, both ways and returns the extended FeatureLocation

        Arguments:
            location: one FeatureLocation instance
            by: integer, to be extended by
            record: parent record of the location

        Returns:
            a new FeatureLocation that will be extended
    """
    # simple features on linear records are straightforward
    if isinstance(location, (FeatureLocation)) and not record.is_circular():
         return FeatureLocation(max(0, location.start - by),
                                min(location.end + by, len(record)))

    bio_start_by = None
    bio_end_by = None
    upstream_part = None
    downstream_part = None
    #
    # todo: # extension would overlap itself on a circular record -> shrink extension to an acceptable value
    # <code here>

    # first part's end (bio_start) on reverse strand of circular record is overflowing the ori -> create an overflow feature
    if location.parts[0].strand == -1 and record.is_circular() and location.parts[0].end + by > len(record):
        bio_start_by = len(record) - location.parts[0].end
        upstream_part = FeatureLocation(0, by - end_by, strand=location.parts[0].strand)
    # first part's start (bio_start) on forward or unspecified strand of a circular record is overflowing the ori -> create an overflow feature
    elif location.parts[0].strand != -1 and record.is_circular() and location.parts[0].start - by < 0:
        bio_start_by = -location.parts[0].start
        upstream_part = FeatureLocation(len(record) - by - bio_start_by, len(record), strand=location.parts[0].strand)
    # else no overflow
    else:
        bio_start_by = by if location.parts[0].strand == -1 else -by

    # last part's start (bio_end) on reverse strand of circular record is overflowing the ori -> create an overflow feature
    if location.parts[-1].strand == -1 and record.is_circular() and location.parts[-1].start - by < 0:
        bio_end_by = -location.parts[-1].start
        downstream_part = FeatureLocation(len(record) - by - bio_end_by, len(record), strand=location.parts[-1].strand)
    # last part's end (bio_end) on forward or unspecified strand of a circular record is overflowing the ori -> create an overflow feature
    if location.parts[-1].strand != -1 and record.is_circular() and location.parts[-1].end + by > len(record):
        bio_end_by = len(record) - location.parts[-1].end
        downstream_part = FeatureLocation(0, by - bio_end_by, strand=location.parts[-1].strand)
    # else no overflow
    else:
        bio_end_by = -by if location.parts[0].strand == -1 else by

    parts = []
    for part in location.parts:
        start = part.start
        if part.start is location.bio_start:
            start = part.start + bio_start_by
        elif part.start is location.bio_end:
            start = part.start + bio_end_by
        end = part.end
        if part.end is location.bio_start:
            end = part.end + bio_start_by
        elif part.end is location.bio_end:
            end = part.end + bio_end_by
        parts.append(FeatureLocation( start, end, strand=part.strand))
    if upstream_part:
        parts.insert(0, upstream_part)
    if downstream_part:
        parts.append(downstream_part)

    return CompoundLocation(parts) if len(parts) > 1 else parts[0]


def merge_within_container(first: FeatureLocation, second: FeatureLocation, container: FeatureLocation) -> FeatureLocation:
    """ Merges two locations, but only if the merge location would still be fully contained in a container location
        Returns:
            a new, merged FeatureLocation, which is similar to FeatureLocation(min(first.start, second.start),max(first.end, second.end)),
            but allows for cross-ori locations, both as first, second or container
    """
    parts = []
    for container_subpart in container.parts: #only one for most cases, two for ori-split
        contained_parts = [part for part in first.parts + second.parts if locations_overlap(container_subpart, part)]
        parts.append(FeatureLocation(min(part.start for part in contained_parts), max(part.end for part in contained_parts)))

    return CompoundLocation(parts) if len(parts) > 1 else parts[0]


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
