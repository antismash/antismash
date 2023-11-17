# License: GNU Affero General Public License v3 or later
# A copy of GNU AGPL v3 should have been included in this software package in LICENSE.txt.

""" Helper functions for location operations """

import logging
from typing import (
    Iterable,
    List,
    Optional,
    Tuple,
    Type,
    TypeVar,
    Union,
)

from Bio.SeqFeature import (
    Position,
    AfterPosition,
    BeforePosition,
    CompoundLocation as _CompoundLocation,
    ExactPosition,
    Location as _Location,
    SimpleLocation as _SimpleLocation,
    SeqFeature,
    UnknownPosition,
)

from .errors import SecmetInvalidInputError


Location = Union["FeatureLocation", "CompoundLocation"]
# a generic for 'B'iopython types while they aren't sharing a common parent class
B = TypeVar("B", _CompoundLocation, _SimpleLocation)
# a generic for the resulting type of the mixin and biopython classes
T = TypeVar("T", bound=Location)


class _LocationMixin(_Location):
    def crosses_origin(self: T, *, allow_reversing: bool = False) -> bool:
        """ Determines if the location would cross the origin of a record.

            Arguments:
                allow_reversing: if True, checks both possible orderings of exons
                                 regardless of strand

            Returns:
                False if the location does not bridge the origin or if the location
                is of indeterminate strand, otherwise True
        """
        return location_bridges_origin(self, allow_reversing=allow_reversing)

    def clone(self: T) -> T:
        """ Clones the location, ensuring that no mutable attributes will be shared
            with the original.

            Returns:
                a new instance
        """
        raise NotImplementedError()

    def clone_with_frameshift(self: T, start: Union[str, int], *, undo: bool = False) -> T:
        """ Generates a new location to represent a frameshift of an existing location.
            Forward strand locations will have their start coordinate lowered.
            Reverse strand locations will have their end coordinate raised.

            Arguments:
                start: a 1-indexed integer or string as per the genbank "codon_start" qualifier
                undo: whether to treat the frameshift as undoing a previous frameshift

            Returns:
                a new location instance of the same type
        """
        return frameshift_location_by_qualifier(self, start, undo=undo)

    def clone_with_offset(self: T, offset: int, *, wrap_point: int = None) -> T:
        """ Creates a new location at the given offset to the original.
            Will not loop over the origin and offsets cannot make locations negative unless the
            origin/wrapping point is provided.

            Arguments:
                offset: the amount to offset
                wrap_point: the origin/coordinate at which locations wrap around

            Returns:
                a new location instance
        """
        return offset_location(self, offset, wrap_point=wrap_point)

    def contains_overlapping_exons(self: T) -> bool:
        """ Returns True if the location contains multiple exons sharing the same stop codon """
        return location_contains_overlapping_exons(self)

    def contains(self: T, other: T) -> bool:
        """ Returns True if this location contains the given location """
        return location_contains_other(self, other)

    def convert_protein_position_to_dna(self: T, start: int, end: int) -> tuple[int, int]:
        """ Convert a protein position to a nucleotide sequence position for use in generating
            new FeatureLocations from existing FeatureLocations and/or CompoundLocations.

            Arguments:
                position: the position in question, must be contained by the location
                location: the location of the related feature, for handling introns/split locations

            Returns:
                an int representing the calculated DNA location
        """
        return convert_protein_position_to_dna(start, end, self)

    def get_distance_to(self: T, other: T, wrap_point: int = None) -> int:
        """ Finds the shortest distance between the two given features, crossing
            the origin if provided.

            Overlapping features are considered to have zero distance.

            Arguments:
                other: the location to get a distance to
                wrap_point: the point at which locations can wrap, if given

            Returns:
                the distance between the two locations
        """
        return get_distance_between_locations(self, other, wrap_point)


class FeatureLocation(_LocationMixin, _SimpleLocation):
    """ A wrapper of biopython's SimpleLocation (previously FeatureLocation) to add extra
        functionality.
    """
    def clone(self: T) -> T:
        return FeatureLocation(self.start, self.end, self.strand)

    @classmethod
    def from_biopython(cls: Type["FeatureLocation"], bio: _SimpleLocation) -> "FeatureLocation":
        """ Constructs an instance from the given biopython FeatureLocation.

            Arguments:
                bio: the biopython location to convert
        """
        return cls(bio.start, bio.end, bio.strand)

    def __contains__(self, value: Union[int, T]) -> bool:
        if isinstance(value, int):
            return super().__contains__(value)
        return self.contains(value)


SimpleLocation = FeatureLocation  # for name mapping purposes between older and newer biopython styles


class CompoundLocation(_LocationMixin, _CompoundLocation):
    """ A wrapper of biopython's CompoundLocation to add extra functionality.
    """
    def clone(self: T) -> T:
        return CompoundLocation(self.parts.copy(), operator=self.operator)

    @classmethod
    def from_biopython(cls: Type["CompoundLocation"], bio: _CompoundLocation) -> "CompoundLocation":
        """ Constructs an instance from the given biopython CompoundLocation.

            Arguments:
                bio: the biopython location to convert
        """
        return cls(bio.parts, operator=bio.operator)

    def __contains__(self, value: Union[int, T]) -> bool:
        if isinstance(value, int):
            return super().__contains__(value)
        return self.contains(value)


def _reduce_parts_to_location(parts: list[FeatureLocation], wrap_point: Optional[int]) -> Location:
    """ Reduces multiple FeatureLocations into a minimal location that may or may
        not cross the origin
    """
    # if it's already reduced as much as it can be, return early
    if len(parts) == 1:
        return parts[0]
    # if the locations cross the origin, then the resulting reduction has to
    # also cross the origin
    temp = CompoundLocation(parts)
    if location_bridges_origin(temp):
        if wrap_point is None:
            raise ValueError("Cannot merge cross-origin location without a wrap point")
        assert wrap_point > 0  # tested before this function is called, but just in case
        lower, upper = split_origin_bridging_location(temp)
        return CompoundLocation([
            FeatureLocation(min(up.start for up in upper), wrap_point, 1),
            FeatureLocation(0, max(low.end for low in lower), 1),
        ])
    # only the simple case remains, where it's not cross origin but still has multiple exons
    temp = CompoundLocation(parts)
    return FeatureLocation(temp.start, temp.end, strand=temp.strand)


def _merge_over_origin(locations: list[Location], wrap_point: int) -> list[Location]:
    """ Merges, where possible, locations that individually don't cross the origin,
        but have a shorter distance over the origin
    """
    new_locations: list[Location] = []
    merged: set[str] = set()
    # first, compact those parts that can be compacted  # TODO, not sure about split with some configs
    upper, lower = _split_sections_around_origin(locations, wrap_point)
    if lower and upper:
        locations = [connect_locations(upper), connect_locations(lower)]
    elif lower:
        locations = [connect_locations(lower)]
    elif upper:
        locations = [connect_locations(upper)]
    # this loop would normally be a do-while loop, but python doesn't support that
    # so use a boolean to require at least one run
    has_run = False
    while len(locations) > 1 and (merged or not has_run):
        has_run = True
        merged.clear()
        new_locations.clear()
        for i, location in enumerate(locations):
            assert len(location.parts) == 1, location
            for other in locations[i:]:
                assert len(other.parts) == 1
                # if it was already merged into a new location, it'll cause issues this iteration
                if str(other) in merged:
                    continue
                # if over the origin is shorter, merge the two into a cross-origin location
                over_origin = get_distance_between_locations(location, other, wrap_point=wrap_point)
                standard = get_distance_between_locations(location, other)
                if over_origin < standard:
                    if other.start < location.start:
                        upper = FeatureLocation(location.start, wrap_point, 1)
                        lower = FeatureLocation(0, other.end, 1)
                    else:
                        upper = FeatureLocation(other.start, wrap_point, 1)
                        lower = FeatureLocation(0, location.end, 1)
                    location = CompoundLocation([upper, lower])
                    merged.add(str(other))
                    break
            if str(location) not in merged:
                new_locations.append(location)
        assert len(new_locations) <= len(locations)
        if not merged:
            break
        assert len(new_locations) == len(locations) - len(merged), f"{new_locations=}\n{locations=}\n{merged=}"
        locations = list(new_locations)
    return locations


def _split_sections_around_origin(locations: list[Location], origin: int,
                                  ) -> tuple[list[Location], list[Location]]:
    """ Separates a list of locations into pre- and post-origin portions """
    pre_chunks = []
    post_chunks = []
    for location in locations:
        if location_bridges_origin(location):
            lower, upper = split_origin_bridging_location(location)
            pre_chunks.append(_reduce_parts_to_location(upper, origin))
            post_chunks.append(_reduce_parts_to_location(lower, origin))
        else:
            location = FeatureLocation(location.start, location.end, 1)
            # check the distance from each side of the location to each side of the record
            if location.start < origin - location.end:
                post_chunks.append(location)
            else:
                pre_chunks.append(location)
    return pre_chunks, post_chunks


def connect_locations(locations: list[Location], wrap_point: int = None) -> Location:
    """ Creates as small a location as possible that fully covers the given locations.
        With a circular record, connections will cross the origin if smaller.

        The resulting location will always be on the forward strand.

        Arguments:
            locations: the locations to connect
            wrap_point: the record length for circular records, otherwise None

        Returns:
            a location, possibly with two parts if it crosses the origin in a circular record
    """
    if not locations:
        raise ValueError("At least one location is required")
    any_cross_origin = any(location_bridges_origin(loc) for loc in locations)

    # linear records can't have cross-origin features at all, so don't try to handle it
    if any_cross_origin and wrap_point is None:
        raise ValueError("Connecting origin-bridging locations requires the record length")

    locations = [_reduce_parts_to_location(location.parts, wrap_point) for location in locations]

    # handle the simplest case first, non-circular inputs
    if wrap_point is None:
        start = min(loc.start for loc in locations)
        end = max(loc.end for loc in locations)
        strand = locations[0].strand if all(loc.strand == locations[0].strand for loc in locations) else None
        return FeatureLocation(start, end, strand=strand)

    assert wrap_point > 0

    # a little more setup is required in the case of circularity with no origin-crossing features,
    # if it would be shorter to take the path over the origin between two locations
    if not any_cross_origin:
        locations = _merge_over_origin(locations, wrap_point=wrap_point)

    # now that merging has happened, if there's only one location left, don't continue
    # otherwise the recursion later on will be infinite
    if len(locations) == 1:
        return locations[0]

    pre_chunks, post_chunks = _split_sections_around_origin(locations, wrap_point)

    # then combine each chunk into a part
    if not pre_chunks:
        if post_chunks == locations:
            wrap_point = None  # don't try the same process all over again with no change
        result = connect_locations(post_chunks)
    elif not post_chunks:
        if pre_chunks == locations:
            wrap_point = None  # don't try the same process all over again with no change
        result = connect_locations(pre_chunks)
    else:
        pre = connect_locations(pre_chunks, wrap_point)
        post = connect_locations(post_chunks, wrap_point)
        assert isinstance(pre, FeatureLocation), pre
        assert isinstance(post, FeatureLocation), post
        # don't build compound locations if one section is fully contained in the other
        if location_contains_other(pre, post) or location_contains_other(post, pre):
            return pre if len(pre) > len(post) else post
        if locations_overlap(pre, post):
            result = connect_locations([pre, post])
            assert result.strand == 1
        else:
            # all created cross-origin locations are forward strand for clarity
            pre.strand = 1
            post.strand = 1
            result = CompoundLocation([pre, post])
            assert result.strand == 1
    return result


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


def build_location_from_others(locations: list[Location]) -> Location:
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


def location_from_biopython(bio: B) -> Location:
    """ Converts the given biopython location instance into a wrapped version
        for more utility.

        Arguments:
            bio: the biopython location to convert

        Returns:
            an instance of Location, matching the type of the location provided
    """
    if isinstance(bio, _CompoundLocation):
        parts = [location_from_biopython(part) for part in bio.parts]
        return CompoundLocation(parts, bio.operator)
    return FeatureLocation(bio.start, bio.end, bio.strand)


def _is_valid_split(lower: list[FeatureLocation], upper: list[FeatureLocation], strand: int) -> bool:
    """ Returns True if the results of a split are valid:
        - mutually exclusive areas covered
        - each section must be ordered correctly for the strand
    """
    if not lower or not upper:
        return False

    # check that both sections cover a mutually exclusive area
    if locations_overlap(connect_locations(lower), connect_locations(upper)):
        return False

    # check that all components in each section are correctly ordered
    for section in [upper, lower]:
        starts = [part.start for part in section]
        if sorted(starts, reverse=(strand == -1)) != starts:
            return False
    return True


def split_origin_bridging_location(location: Location) -> tuple[
                                                      List[FeatureLocation], List[FeatureLocation]]:
    """ Splits a CompoundLocation into two sections.
        The first contains the low-position parts (immediately after the origin
        in a forward direction), the second handles the high-position parts.

        Arguments:
            location: the CompoundLocation to split

        Returns:
            a tuple of lists, each list containing one or more FeatureLocations
    """
    if isinstance(location, FeatureLocation):
        return ([location], [])

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


def make_forwards(location: Location) -> Location:
    """ Creates a copy of a location in the forward strand,
        reordering the components if it was on the reverse strand.

        Arguments:
            location: the location to convert

        Returns:
            a new location in the forward strand
    """
    parts = [FeatureLocation(part.start, part.end, 1) for part in location.parts]
    if location.strand == -1:
        parts.reverse()
    if len(parts) == 1:
        return parts[0]
    loc = CompoundLocation(parts)
    return loc


def location_contains_overlapping_exons(location: Union[Location, B]) -> bool:
    """ Checks for multiple exons with the same end location, meaning they use the
        same stop codon

        Arguments:
            location: the location to check

        Returns:
            True if the location contains exons sharing a stop codon
    """
    if not isinstance(location, (CompoundLocation, FeatureLocation, _CompoundLocation, _SimpleLocation)):
        raise TypeError(f"expected location type, received {type(location)}")
    if len(location.parts) == 1:
        return False
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
            raise ValueError(f"feature outside record sequence: {feature.location}")
        # features with overlapping exons cause translation problems
        if location_contains_overlapping_exons(feature.location):
            raise ValueError(f"location contains overlapping exons: {feature.location}")

    # non-circular records with compound locations need to have the right part ordering
    # for translations, only really relevant for reverse strand features
    # first find what pattern has been used for locations
    standard = 0
    non_standard = 0
    for feature in features:
        # update from biopython to internal types
        feature.location = location_from_biopython(feature.location)

        if not feature.location.strand or feature.type not in ["CDS", "gene"]:
            continue

        if location_bridges_origin(feature.location):
            non_standard += 1
        else:
            standard += 1

    if can_be_circular:
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


def offset_location(location: Location, offset: int, *, wrap_point: int = None) -> Location:
    """ Creates a new location at the given offset to the original.
        Will not loop over the origin and offsets cannot make locations negative unless the
        origin/wrapping point is provided.

        Arguments:
            location: the location to shift
            offset: the amount to offset
            wrap_point: the origin/coordinate at which locations wrap around

        Returns:
            a new location instance
    """
    def shifted_location() -> Location:
        if not offset:
            return location.clone()
        parts = location.parts
        new = []
        for part in parts:
            start = part.start + offset
            end = part.end + offset
            assert start < end
            assert wrap_point is not None or start >= 0 and end > 0
            new.append(FeatureLocation(start, end, strand=part.strand))
        if isinstance(location, CompoundLocation):
            return CompoundLocation(new, operator=location.operator)
        return new[0]

    # the trivial cases can be handled easily
    if not wrap_point or not offset:
        return shifted_location()

    if wrap_point < 1:
        raise ValueError(f"wrapping point must be positive: {wrap_point}")

    # if the location covered the entire area, don't adjust it at all, since it'll be the same
    if len(location) == wrap_point:
        return location.clone()

    # the trivial case, no wrapping required
    if 0 < location.start + offset < location.end + offset < wrap_point:
        return shifted_location()

    # the remaining cases may either start over the origin or result in being over the origin
    # start by just shifting everything along, splitting where required
    parts = shifted_location().parts

    new_parts = []
    for part in parts:
        # the shifted locations may be over the record or before the record,
        # so adjust as necessary
        start = (part.start + wrap_point) % wrap_point
        end = (part.end + wrap_point) % wrap_point
        # if the current part is still within the record, it can be used as is
        if 0 <= start < end <= wrap_point:
            new_parts.append(FeatureLocation(start, end, part.strand))
            continue
        # otherwise split the shifted part if it now crosses the origin
        new_parts.extend([
            FeatureLocation(start, wrap_point, part.strand),
            FeatureLocation(0, end, part.strand),
        ])
    for part in new_parts:
        assert 0 <= part.start < part.end <= wrap_point, part

    # if originally cross-origin, then a merge may be required
    previous = new_parts[0]
    merged = [previous]
    for part in new_parts[1:]:
        if previous.end == part.start:
            assert previous.strand == part.strand
            # replace the existing one
            merged[-1] = FeatureLocation(previous.start, part.end, part.strand)
        else:
            merged.append(part)
        previous = part

    assert merged

    return CompoundLocation(merged) if len(merged) > 1 else merged[0]


def remove_redundant_exons(location: Location) -> Location:
    """ Generates a new location that with redudant exons removed.
        Redundant exons are those that cover a location already covered by a larger exon.

        Arguments:
            location: the location to trim

        Returns:
            a new location instance, if redundant exons are found, otherwise the existing location
    """
    if isinstance(location, FeatureLocation):
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
