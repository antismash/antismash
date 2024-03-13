# License: GNU Affero General Public License v3 or later
# A copy of GNU AGPL v3 should have been included in this software package in LICENSE.txt.

# for test files, silence irrelevant and noisy pylint warnings
# pylint: disable=use-implicit-booleaness-not-comparison,protected-access,missing-docstring
# reusing variables in different argument positions seems to trick pylint
# pylint: disable=arguments-out-of-order

import unittest
from unittest.mock import patch

from antismash.common import secmet
from antismash.common.secmet.locations import (
    _adjust_location_by_offset as adjust,
    convert_protein_position_to_dna,
    build_location_from_others,
    connect_locations,
    ensure_valid_locations,
    get_distance_between_locations as _get_distance_between_locations,
    location_bridges_origin,
    split_origin_bridging_location as splitter,
    location_contains_other as _location_contains_other,
    location_contains_overlapping_exons,
    location_from_string,
    locations_overlap,
    offset_location as _offset_location,
    remove_redundant_exons,
    FeatureLocation,
    _SimpleLocation as BioSimple,
    _CompoundLocation as BioCompound,
    CompoundLocation,
    AfterPosition,
    BeforePosition,
    ExactPosition,
    SeqFeature,
    UnknownPosition,
)


class TestConnectLocations(unittest.TestCase):
    def setUp(self):
        self.func = connect_locations

    def test_finding_shortest_path(self):
        locations = [FeatureLocation(5, 10, 1), FeatureLocation(90, 95, 1)]
        # with the origin not set, the record isn't circular and the joined location
        # will be the long way around, covering most of the record
        expected_simple = FeatureLocation(5, 95, 1)
        result = self.func(locations)
        assert result == expected_simple
        # but when the origin is set, and would provide a shorter path, the joined
        # location will cross the origin
        result = self.func(locations, wrap_point=100)
        assert result == CompoundLocation([FeatureLocation(90, 100, 1), FeatureLocation(0, 10, 1)])
        # and to check that the distance calculation actually chooses the shortest path
        result = self.func(locations, wrap_point=200)
        assert result == expected_simple

    def test_cross_origin_inputs(self):
        # when two locations cross the origin, they should still end up with just two parts
        expected = CompoundLocation([FeatureLocation(7644228, 7663439, 1), FeatureLocation(0, 9008, 1)])
        result = self.func([
            CompoundLocation([FeatureLocation(7644228, 7663439, 1), FeatureLocation(0, 19, 1)]),
            CompoundLocation([FeatureLocation(7663426, 7663439, 1), FeatureLocation(0, 9008, 1)]),
        ], wrap_point=7663439)
        assert result == expected

    def test_cross_origin_inputs_reverse(self):
        # when two locations cross the origin, they should still end up with just two parts
        expected = CompoundLocation([FeatureLocation(7644228, 7663439, 1), FeatureLocation(0, 9008, 1)])
        result = self.func([
            CompoundLocation([FeatureLocation(0, 19, -1), FeatureLocation(7644228, 7663439, -1)]),
            CompoundLocation([FeatureLocation(0, 9008, -1), FeatureLocation(7663426, 7663439, -1)]),
        ], wrap_point=7663439)
        assert result == expected

    def test_cross_origin_in_linear(self):
        location = CompoundLocation([FeatureLocation(760, 766, 1), FeatureLocation(0, 90, 1)])
        with self.assertRaisesRegex(ValueError, "origin-bridging .* requires the record length"):
            self.func([location])

    def test_contained(self):
        # if one of the ocations is contained by a cross-origin location, it shouldn't break
        expected = CompoundLocation([FeatureLocation(760, 766, 1), FeatureLocation(0, 90, 1)])
        result = self.func([expected, expected.parts[1]], wrap_point=766)
        assert result == expected

    def test_most_of_area(self):
        record_length = 98206
        locations = [
            CompoundLocation([
                FeatureLocation(65398, record_length, 1),
                FeatureLocation(0, 52159, 1),
            ]),
            CompoundLocation([
                FeatureLocation(93497, record_length, 1),
                FeatureLocation(0, 52159, 1),
            ]),
            FeatureLocation(7053, 52159, 1),  # more than half of the record
        ]
        result = self.func(locations, wrap_point=record_length)
        assert result.crosses_origin()
        assert result == CompoundLocation([
            FeatureLocation(65398, record_length, 1),
            FeatureLocation(0, 52159, 1),
        ])

    def test_overlap(self):
        locations = [FeatureLocation(0, 60, 1), FeatureLocation(50, 100, 1)]
        # both with and without circularity, the end result should be a non-compound location
        # covering the full record
        expected = FeatureLocation(0, 100, 1)
        for record_length in [10, 100]:
            result = self.func(locations, wrap_point=record_length)
            assert result == expected
        assert self.func(locations) == expected
        with self.assertRaises(AssertionError):
            self.func(locations, wrap_point=0)

    def test_overlaps(self):
        locations = [
            CompoundLocation([
                FeatureLocation(85_241, 172_542, strand=1),
                FeatureLocation(0, 85_240, strand=1)
            ], 'join'),
            CompoundLocation([
                FeatureLocation(89_576, 172_542, strand=1),
                FeatureLocation(0, 89_575, strand=1)
            ], 'join'),
            FeatureLocation(47_775, 96_472, strand=1),
            FeatureLocation(57_556, 100_658, strand=1),
            FeatureLocation(78_493, 125_885, strand=1)
        ]
        wrap_point = 172_542
        result = self.func(locations, wrap_point=wrap_point)
        assert result == FeatureLocation(0, 172_542, strand=1)

    def test_odd(self):
        locations = [
            FeatureLocation(331_108, 332_344, strand=1),
            FeatureLocation(24, 831, strand=1),
            FeatureLocation(975, 1137, strand=1),
            FeatureLocation(1764, 2247, strand=1),
        ]
        wrap_point = 340_000
        result = self.func(locations, wrap_point=wrap_point)
        assert len(result.parts) == 2
        assert result.parts[0] == FeatureLocation(331_108, 340_000, 1)
        assert result.parts[1] == FeatureLocation(0, 2247, 1)

class TestProteinPositionConversion(unittest.TestCase):
    def func(self, start, end, location):
        static = convert_protein_position_to_dna(start, end, location)
        dynamic = location.convert_protein_position_to_dna(start, end)
        assert static == dynamic
        return static

    def test_position_conversion_simple_forward(self):
        location = FeatureLocation(0, 15, strand=1)
        assert len(location) == 15
        assert self.func(0, 2, location) == (0, 6)
        assert self.func(1, 4, location) == (3, 12)

    def test_position_conversion_simple_reverse(self):
        location = FeatureLocation(0, 15, strand=-1)
        assert len(location) == 15
        assert self.func(0, 2, location) == (9, 15)
        assert self.func(1, 4, location) == (3, 12)

    def test_position_conversion_nonzero_start(self):
        location = FeatureLocation(6, 21, strand=1)
        assert len(location) == 15
        assert self.func(0, 2, location) == (6, 12)
        assert self.func(1, 4, location) == (9, 18)

        location = FeatureLocation(6, 21, strand=-1)
        assert len(location) == 15
        assert self.func(0, 2, location) == (15, 21)
        assert self.func(1, 4, location) == (9, 18)

    def test_position_conversion_nonzero_compound(self):
        location = CompoundLocation([FeatureLocation(6, 18, strand=1),
                                     FeatureLocation(24, 27, strand=1)])
        assert len(location) == 15
        assert self.func(0, 2, location) == (6, 12)
        assert self.func(1, 4, location) == (9, 18)
        assert self.func(3, 5, location) == (15, 27)

        location = CompoundLocation([FeatureLocation(6, 15, strand=-1),
                                     FeatureLocation(21, 27, strand=-1)])
        assert len(location) == 15
        assert self.func(0, 2, location) == (21, 27)
        assert self.func(1, 4, location) == (9, 24)
        assert self.func(3, 5, location) == (6, 12)

    def test_position_conversion_compound_forward(self):
        location = CompoundLocation([FeatureLocation(0, 6, strand=1),
                                     FeatureLocation(9, 18, strand=1)])
        assert len(location) == 15
        assert self.func(0, 4, location) == (0, 15)
        assert self.func(1, 5, location) == (3, 18)

        location = CompoundLocation([FeatureLocation(0, 6, strand=1),
                                     FeatureLocation(12, 15, strand=1),
                                     FeatureLocation(21, 27, strand=1)])
        assert len(location) == 15
        assert self.func(0, 4, location) == (0, 24)
        assert self.func(1, 5, location) == (3, 27)
        assert self.func(2, 3, location) == (12, 15)

    def test_position_conversion_compound_reverse(self):
        location = CompoundLocation([FeatureLocation(0, 6, strand=-1),
                                     FeatureLocation(9, 18, strand=-1)])
        assert len(location) == 15
        assert self.func(0, 4, location) == (3, 18)
        assert self.func(1, 5, location) == (0, 15)

        location = CompoundLocation([FeatureLocation(0, 6, strand=-1),
                                     FeatureLocation(12, 15, strand=-1),
                                     FeatureLocation(21, 27, strand=-1)])
        assert len(location) == 15
        assert self.func(0, 4, location) == (3, 27)
        assert self.func(1, 5, location) == (0, 24)
        assert self.func(2, 3, location) == (12, 15)

    def test_other(self):
        location = CompoundLocation([FeatureLocation(5922, 6190, strand=1),
                                     FeatureLocation(5741, 5877, strand=1),
                                     FeatureLocation(4952, 5682, strand=1)])
        assert self.func(97, 336, location) == (5243, 6064)

        location = CompoundLocation([FeatureLocation(5922, 6190, strand=-1),
                                     FeatureLocation(5741, 5877, strand=-1),
                                     FeatureLocation(4952, 5682, strand=-1)])
        assert self.func(97, 336, location) == (5078, 5854)


class TestCompoundCombination(unittest.TestCase):
    def test_separate(self):
        locations = [FeatureLocation(6, 9, 1), FeatureLocation(12, 16, 1)]
        new = build_location_from_others(locations)
        assert isinstance(new, CompoundLocation)
        assert new == CompoundLocation(locations)

    def test_single(self):
        for strand in [1, -1]:
            location = FeatureLocation(6, 9, strand)
            new = build_location_from_others([location])
            assert isinstance(new, FeatureLocation) and not isinstance(new, CompoundLocation)
            assert new == location

    def test_single_compound(self):
        for strand in [1, -1]:
            location = CompoundLocation([FeatureLocation(6, 9, strand),
                                         FeatureLocation(12, 16, strand)])
            new = build_location_from_others([location])
            assert new == location

    def test_all_merged(self):
        for strand in [1, -1]:
            locations = [FeatureLocation(6, 9, strand),
                         FeatureLocation(9, 12, strand),
                         FeatureLocation(12, 16, strand)]
            new = build_location_from_others(locations)
            assert isinstance(new, FeatureLocation) and not isinstance(new, CompoundLocation)
            assert new == FeatureLocation(6, 16, strand)

    def test_some_merged(self):
        for strand in [1, -1]:
            locations = [FeatureLocation(1, 4, strand),
                         FeatureLocation(6, 9, strand),
                         FeatureLocation(9, 12, strand),
                         FeatureLocation(15, 18, strand)]
            new = build_location_from_others(locations)
            assert isinstance(new, CompoundLocation)
            assert new == CompoundLocation([FeatureLocation(1, 4, strand),
                                            FeatureLocation(6, 12, strand),
                                            FeatureLocation(15, 18, strand)])


def build_compound(pairs, strand, operator="join"):
    assert len(pairs) >= 2, "invalid CompoundLocation would be created"
    parts = []
    for start, end in pairs:
        parts.append(FeatureLocation(start, end, strand))
    return CompoundLocation(parts, operator=operator)


def is_bridged(location, **kwargs):
    static = location_bridges_origin(location, **kwargs)
    dynamic = location.crosses_origin(**kwargs)
    assert static == dynamic
    return static


class TestBridgeDetection(unittest.TestCase):
    def test_forward(self):
        assert is_bridged(build_compound([(9, 12), (0, 3)], 1))
        assert is_bridged(build_compound([(9, 12), (0, 3), (4, 5)], 1))
        assert is_bridged(build_compound([(4, 5), (9, 12), (0, 3)], 1))
        assert not is_bridged(build_compound([(0, 3), (9, 12)], 1))

    def test_reverse(self):
        assert is_bridged(build_compound([(0, 3), (9, 12)], -1))
        assert is_bridged(build_compound([(6, 9), (0, 3), (15, 18)], -1))
        assert is_bridged(build_compound([(0, 3), (15, 18), (6, 9)], -1))
        assert not is_bridged(build_compound([(9, 12), (0, 3)], -1))

    def test_alternate_orderings(self):
        loc = build_compound([(0, 3), (6, 9), (12, 15)], -1)
        assert loc.parts[0].start == 0
        assert is_bridged(loc)
        assert loc.parts[0].start == 0

        assert not is_bridged(loc, allow_reversing=True)
        assert loc.parts[0].start == 12

        loc = build_compound([(12, 15), (6, 9), (0, 3)], -1)
        assert loc.parts[0].start == 12
        assert not is_bridged(loc, allow_reversing=True)
        assert loc.parts[0].start == 12

    def test_bad_strand(self):
        pairs = [(9, 12), (0, 3)]
        assert is_bridged(build_compound(pairs, 1))
        assert not is_bridged(build_compound(pairs, None))

    def test_not_bridged(self):
        assert not is_bridged(build_compound([(1, 6), (5, 10)], 1))
        assert not is_bridged(build_compound([(5, 10), (1, 6)], -1))

    def test_indeterminate(self):
        assert is_bridged(build_compound([(0, 3), (12, 15), (6, 9)], -1), allow_reversing=True)


class TestBridgedSplit(unittest.TestCase):
    def check_pairs(self, parts, pairs):
        assert [(int(part.start), int(part.end)) for part in parts] == pairs

    def test_simple_forward(self):
        loc = build_compound([(9, 12), (0, 3)], 1)
        lower, upper = splitter(loc)
        self.check_pairs(lower, [(0, 3)])
        self.check_pairs(upper, [(9, 12)])

    def test_simple_reverse(self):
        loc = build_compound([(0, 3), (9, 12)], -1)
        lower, upper = splitter(loc)
        self.check_pairs(lower, [(0, 3)])
        self.check_pairs(upper, [(9, 12)])

    def test_extras_forward(self):
        loc = build_compound([(15, 18), (0, 3), (6, 9)], 1)
        lower, upper = splitter(loc)
        self.check_pairs(lower, [(0, 3), (6, 9)])
        self.check_pairs(upper, [(15, 18)])

        loc = build_compound([(6, 9), (15, 18), (0, 3)], 1)
        lower, upper = splitter(loc)
        self.check_pairs(lower, [(0, 3)])
        self.check_pairs(upper, [(6, 9), (15, 18)])

    def test_extras_reverse(self):
        loc = build_compound([(6, 9), (0, 3), (15, 18)], -1)
        lower, upper = splitter(loc)
        self.check_pairs(lower, [(6, 9), (0, 3)])
        self.check_pairs(upper, [(15, 18)])

        loc = build_compound([(0, 3), (15, 18), (6, 9)], -1)
        lower, upper = splitter(loc)
        self.check_pairs(lower, [(0, 3)])
        self.check_pairs(upper, [(15, 18), (6, 9)])

    def test_unusable(self):
        # this format crosses the origin multiple times and can't be interpreted
        raw_loc_parts = [(0, 3), (6, 9), (12, 15), (18, 21)]
        # test a variety of badly ordered parts
        for ordering in [[0, 2, 1, 3], [0, 1, 3, 2], [0, 3, 1, 2]]:
            # cycle them around to test position independence
            for i in range(len(ordering)):
                loc_parts = [raw_loc_parts[i] for i in ordering[i:] + ordering[:i]]

                # forward
                loc = build_compound(loc_parts, 1)
                assert is_bridged(loc)
                with self.assertRaisesRegex(ValueError, "cannot determine correct ordering"):
                    splitter(loc)

                # reverse
                loc = build_compound(loc_parts[::-1], -1)
                assert is_bridged(loc)
                with self.assertRaisesRegex(ValueError, "cannot determine correct ordering"):
                    splitter(loc)

    def test_not_bridging_forward(self):
        loc = build_compound([(0, 3), (9, 12)], 1)
        with self.assertRaisesRegex(ValueError, "Location does not bridge origin"):
            print(splitter(loc))

    def test_not_bridging_reverse(self):
        loc = build_compound([(9, 12), (0, 3)], -1)
        with self.assertRaisesRegex(ValueError, "Location does not bridge origin"):
            print(splitter(loc))

        loc = build_compound([(15, 18), (9, 12), (0, 3)], -1)
        with self.assertRaisesRegex(ValueError, "Location does not bridge origin"):
            print(splitter(loc))

    def test_bad_strand(self):
        loc = build_compound([(9, 12), (0, 3)], -1)
        loc.parts[0].strand = 1
        loc.parts[1].strand = -1
        assert loc.strand is None
        with self.assertRaisesRegex(ValueError, "Cannot separate bridged location without a valid strand"):
            print(splitter(loc))

    def test_frameshifts(self):
        loc = build_compound([(4772224, 4772573), (4772572, 4773186), (0, 258)], 1)
        low, high = splitter(loc)
        assert low == loc.parts[2:]
        assert high == loc.parts[0:2]


class TestLocationSerialiser(unittest.TestCase):
    def convert(self, location, expected_type=FeatureLocation):
        assert isinstance(location, expected_type)

        before_string = str(location)
        print(before_string)  # just for help when debugging a failing test
        after_string = str(location)
        assert isinstance(after_string, str)
        assert before_string == after_string

        new_location = location_from_string(after_string)
        assert isinstance(new_location, expected_type)

        return new_location

    def test_before_position(self):
        location = FeatureLocation(BeforePosition(1), ExactPosition(6), strand=-1)
        new_location = self.convert(location)

        assert isinstance(new_location.start, BeforePosition)
        assert new_location.start == 1

        assert isinstance(new_location.end, ExactPosition)
        assert new_location.end == 6

    def test_after_position(self):
        location = FeatureLocation(ExactPosition(1), AfterPosition(6), strand=1)
        new_location = self.convert(location)

        assert isinstance(new_location.start, ExactPosition)
        assert new_location.start == 1

        assert isinstance(new_location.end, AfterPosition)
        assert new_location.end == 6

    def test_unknown_position(self):
        location = FeatureLocation(ExactPosition(1), UnknownPosition(), strand=1)
        new_location = self.convert(location)

        assert isinstance(new_location.start, ExactPosition)
        assert new_location.start == 1

        assert isinstance(new_location.end, UnknownPosition)

    def test_compound(self):
        first = FeatureLocation(1, 6, strand=1)
        second = FeatureLocation(10, 16, strand=1)
        location = CompoundLocation([first, second], operator="join")
        assert 5 in location
        assert 7 not in location
        assert 15 in location

        new_location = self.convert(location, expected_type=CompoundLocation)
        assert location.start == 1
        assert 5 in new_location
        assert 7 not in new_location
        assert 15 in new_location
        assert location.end == 16
        assert new_location.operator == "join"

    def test_strands(self):
        for strand in [1, 0, -1, None]:
            location = FeatureLocation(1, 6, strand=strand)
            new_location = self.convert(location)
            assert new_location.strand == strand


class TestOverlaps(unittest.TestCase):
    def test_simple_simple(self):
        assert not locations_overlap(FeatureLocation(1, 5, strand=1), FeatureLocation(10, 15, strand=1))
        assert locations_overlap(FeatureLocation(1, 25, strand=1), FeatureLocation(10, 15, strand=1))
        assert locations_overlap(FeatureLocation(1, 12, strand=1), FeatureLocation(10, 15, strand=1))

        assert locations_overlap(FeatureLocation(12, 22, strand=-1), FeatureLocation(10, 15, strand=1))
        assert not locations_overlap(FeatureLocation(12, 22, strand=-1), FeatureLocation(10, 12, strand=1))

    def test_mixed(self):
        compound = build_compound([(0, 10), (20, 30), (40, 50)], strand=1)
        simple = FeatureLocation(15, 17)
        assert not locations_overlap(simple, compound)
        assert not locations_overlap(compound, simple)

        simple = FeatureLocation(22, 25)
        assert locations_overlap(simple, compound)
        assert locations_overlap(compound, simple)

        simple = FeatureLocation(35, 45)
        assert locations_overlap(simple, compound)
        assert locations_overlap(compound, simple)

    def test_compound_compound(self):
        first = build_compound([(0, 10), (20, 30), (40, 50)], strand=1)
        second = build_compound([(12, 18), (32, 38), (52, 58)], strand=1)
        assert not locations_overlap(first, second)
        assert not locations_overlap(second, first)

        second = build_compound([(12, 18), (28, 38), (52, 58)], strand=1)
        assert locations_overlap(first, second)
        assert locations_overlap(second, first)

        second = build_compound([(12, 18), (32, 38), (42, 58)], strand=-1)
        assert locations_overlap(first, second)
        assert locations_overlap(second, first)


def location_contains_other(outer, inner):
    static = _location_contains_other(outer, inner)
    dynamic = outer.contains(inner)
    assert isinstance(outer, (CompoundLocation, FeatureLocation))
    assert isinstance(inner, (CompoundLocation, FeatureLocation))
    builtin = inner in outer
    assert static == dynamic == builtin
    return static


class TestContainsOther(unittest.TestCase):
    def test_simple_in_simple(self):
        inner = FeatureLocation(5, 10)
        outer = FeatureLocation(1, 20)

        # clear contains
        assert location_contains_other(outer, inner)
        assert not location_contains_other(inner, outer)

        # on one edge
        outer = FeatureLocation(5, 20)
        assert location_contains_other(outer, inner)
        assert not location_contains_other(inner, outer)

        # on both edges
        outer = FeatureLocation(1, 20)
        assert location_contains_other(outer, inner)
        assert not location_contains_other(inner, outer)

    def test_simple_in_compound(self):
        simple = FeatureLocation(5, 10)
        compound = build_compound([(1, 4), (12, 20)], strand=1)
        assert not location_contains_other(compound, simple)

        simple = FeatureLocation(1, 20)
        assert not location_contains_other(compound, simple)

        simple = FeatureLocation(15, 18)
        assert location_contains_other(compound, simple)

        for part in compound.parts:
            assert location_contains_other(compound, part)

    def test_compound_in_simple(self):
        simple = FeatureLocation(10, 40)
        compound = build_compound([(10, 20), (20, 40)], strand=1)
        assert location_contains_other(simple, compound)

        compound = build_compound([(10, 20), (20, 40), (50, 60)], strand=1)
        assert not location_contains_other(simple, compound)


def overlapping_exons(location):
    static = location_contains_overlapping_exons(location)
    dynamic = location.contains_overlapping_exons()
    assert static == dynamic
    return static


class TestOverlappingExons(unittest.TestCase):
    def test_non_compound(self):
        assert not overlapping_exons(FeatureLocation(10, 40))

    def test_not_overlapping(self):
        assert not overlapping_exons(build_compound([(10, 30), (40, 70)], strand=1))

    def test_overlapping(self):
        assert overlapping_exons(build_compound([(10, 30), (20, 30)], strand=1))
        assert overlapping_exons(build_compound([(10, 30), (70, 100), (20, 30)], strand=1))
        assert overlapping_exons(build_compound([(70, 100), (20, 30), (10, 30)], strand=1))

    def test_bad_types(self):
        for bad in [None, "loc", [FeatureLocation(10, 40)], 5]:
            with self.assertRaises(TypeError):
                overlapping_exons(bad)


class TestEnsureValid(unittest.TestCase):
    def check(self, features, circular=False, seq_len=100):
        return ensure_valid_locations(features, circular, seq_len)

    def test_reversed(self):
        features = [SeqFeature(build_compound([(10, 30), (40, 70)], -1), type="CDS"),
                    SeqFeature(build_compound([(10, 30), (40, 70)], -1), type="gene"),
                    SeqFeature(build_compound([(40, 70), (10, 30)], None), type="other")]
        self.check(features)
        for feature in features:
            assert feature.location.parts[0].start == 40

    def test_inconsistent(self):
        features = [SeqFeature(build_compound([(10, 30), (40, 70)], -1), type="CDS"),
                    SeqFeature(build_compound([(40, 70), (10, 30)], -1), type="gene"),
                    SeqFeature(build_compound([(40, 70), (10, 30)], None), type="other")]
        with self.assertRaisesRegex(ValueError, "inconsistent exon ordering"):
            self.check(features)

    def test_outside_seq(self):
        features = [SeqFeature(FeatureLocation(50, 140, 1))]
        with self.assertRaisesRegex(ValueError, "feature outside record sequence"):
            self.check(features)

    @patch.object(secmet.locations, "location_contains_overlapping_exons", return_value=True)
    def test_overlapping_exons(self, _patched_overlap):
        features = [SeqFeature(FeatureLocation(5, 8, 1))]
        with self.assertRaisesRegex(ValueError, "contains overlapping exons"):
            self.check(features)

    @patch.object(secmet.locations, "location_bridges_origin", return_value=True)
    def test_too_many_in_circular(self, _patched_bridge):
        features = [SeqFeature(build_compound([(10, 30), (0, 9)], -1), type="CDS"),
                    SeqFeature(build_compound([(10, 30), (0, 9)], -1), type="gene"),
                    SeqFeature(build_compound([(10, 30), (0, 9)], -1), type="CDS")]
        with self.assertRaisesRegex(ValueError, "inconsistent exon ordering"):
            self.check(features, circular=True)


class TestLocationAdjustment(unittest.TestCase):
    def setUp(self):
        # not all of these really make sense biologically, but they're all valid computationally
        self.position_types = [ExactPosition, AfterPosition, BeforePosition, int]

    def test_single_forward(self):
        for position_type in self.position_types:
            old = FeatureLocation(position_type(5), 12, 1)
            for offset in range(-2, 3):
                new = adjust(old, offset)
                assert isinstance(new.start, position_type)
                assert new.start == old.start + offset
                assert new.end is old.end

    def test_single_reverse(self):
        for position_type in self.position_types:
            old = FeatureLocation(5, position_type(12), -1)
            for offset in range(-2, 3):
                new = adjust(old, offset)
                assert isinstance(new.end, position_type)
                assert new.end == old.end + offset
                assert new.start is old.start

    def test_compound_forward(self):
        for position_type in self.position_types:
            old = CompoundLocation([FeatureLocation(position_type(5), 12, 1),
                                    FeatureLocation(15, 17, 1)])
            for offset in range(-2, 3):
                new = adjust(old, offset)
                assert isinstance(new.start, position_type)
                assert new.start == old.start + offset
                assert new.parts[0].end is old.parts[0].end
                for old_part, new_part in zip(old.parts[1:], new.parts[1:]):
                    assert old_part is new_part

    def test_compound_reverse(self):
        for position_type in self.position_types:
            old = CompoundLocation([FeatureLocation(15, position_type(17), -1),
                                    FeatureLocation(5, 12, -1)])
            for offset in range(-2, 3):
                new = adjust(old, offset)
                assert isinstance(new.end, position_type)
                assert new.end == old.end + offset
                assert new.parts[0].start is old.parts[0].start
                for old_part, new_part in zip(old.parts[1:], new.parts[1:]):
                    assert old_part is new_part


def offset_location(location, offset, **kwargs):
    static = _offset_location(location, offset, **kwargs)
    dynamic = location.clone_with_offset(offset, **kwargs)
    assert static == dynamic
    return static


class TestOffset(unittest.TestCase):
    def test_simple(self):
        for strand in [-1, None, 1]:
            old = FeatureLocation(5, 10, strand=strand)
            for offset in [-2, 0, 15]:
                new = offset_location(old, offset)
                assert new.start == old.start + offset
                assert new.end == old.end + offset
                assert new.strand == old.strand

    def test_negative(self):
        old = FeatureLocation(5, 10)
        with self.assertRaises(AssertionError):
            offset_location(old, -old.end)

    def test_compound(self):
        for strand in [-1, None, 1]:
            parts = [FeatureLocation(5, 10, strand=strand), FeatureLocation(15, 20, strand=strand)]
            for operator in ["join", "order"]:
                old = CompoundLocation(parts, operator=operator)
                for offset in [-2, 0, 15]:
                    new = offset_location(old, offset)
                    assert len(old.parts) == len(new.parts)
                    assert new.start == old.start + offset
                    assert new.end == old.end + offset
                    assert new.strand == old.strand
                    assert new.operator == old.operator

    def test_wrapping_starting_compound(self):
        loc = CompoundLocation([FeatureLocation(54616, 55016, 1), FeatureLocation(0, 1403, 1)])
        new = offset_location(loc, -41624, wrap_point=55016)
        assert new == FeatureLocation(12992, 14795, 1)

        # offsets in either direction that land on the wrapping point need to be equivalent
        expected = FeatureLocation(0, 200, 1)
        loc = CompoundLocation([FeatureLocation(400, 500, 1), FeatureLocation(0, 100, 1)])
        for offset in [-400, 100]:
            new = offset_location(loc, offset, wrap_point=500)
            assert new == expected

    def test_wrapping_full_width(self):
        full_len = 55016
        loc = FeatureLocation(0, full_len, 1)
        new = offset_location(loc, - full_len // 2, wrap_point=full_len)
        assert new == loc

    def test_wrapping_simple(self):
        loc = FeatureLocation(5, 15, 1)
        assert offset_location(loc, 5, wrap_point=30) == FeatureLocation(10, 20, 1)
        loc = FeatureLocation(5, 15, 1)
        assert offset_location(loc, -20, wrap_point=30) == FeatureLocation(15, 25, 1)

    def test_newly_wrapped(self):
        loc = FeatureLocation(5, 15, 1)
        new = offset_location(loc, 20, wrap_point=30)
        expected = CompoundLocation([
            FeatureLocation(25, 30, 1),
            FeatureLocation(0, 5, 1),
        ])
        assert new == expected

        # and with negative offset
        new = offset_location(loc, -10, wrap_point=20)
        expected = CompoundLocation([
            FeatureLocation(15, 20, 1),
            FeatureLocation(0, 5, 1),
        ])
        assert new == expected

    def test_wrap_point(self):
        loc = FeatureLocation(5, 10, 1)
        for bad in [-5]:
            with self.assertRaisesRegex(ValueError, "must be positive"):
                offset_location(loc, 10, wrap_point=bad)

    def test_offset_also_cross_origin(self):
        # where a location already crosses the origin but the offset doesn't change that
        old = CompoundLocation([
            FeatureLocation(95, 100, 1),
            FeatureLocation(0, 5, 1)
        ])
        new = offset_location(old, 2, wrap_point=old.end)
        assert new.crosses_origin()
        assert new.parts[0].start == 97
        assert new.parts[1].end == 7

    def test_zero_offset(self):
        old = CompoundLocation([
            FeatureLocation(95, 100, 1),
            FeatureLocation(0, 5, 1)
        ])
        new = offset_location(old, 0, wrap_point=old.end)
        assert old == new


class TestRemoveRedundant(unittest.TestCase):
    def create(self, raw_parts, strand=1, operator="join"):
        parts = [FeatureLocation(part[0], part[1], strand=strand) for part in raw_parts]
        if len(parts) == 1:
            return parts[0]
        return CompoundLocation(parts, operator)

    def check_ordering(self, strand=1, operator="join"):
        parts = [(6, 9), (15, 21), (3, 12)]
        for i in range(len(parts)):
            old = self.create(parts[i:] + parts[:i], strand=strand, operator=operator)
            new = remove_redundant_exons(old)
            assert new != old
            assert len(new.parts) == 2
            assert new.parts == [part for part in old.parts if part.start != 6]

    def test_non_compound(self):
        old = self.create([(3, 6)])
        new = remove_redundant_exons(old)
        assert new == old

    def test_non_redundant(self):
        # completely separate
        old = self.create([(3, 6), (9, 12)])
        new = remove_redundant_exons(old)
        assert new == old

        # incomplete overlap
        old = self.create([(3, 6), (4, 7)])
        new = remove_redundant_exons(old)
        assert new == old

    def test_compound_to_non_compound(self):
        old = self.create([(3, 12), (6, 9)])
        new = remove_redundant_exons(old)
        assert not isinstance(new, CompoundLocation)
        assert old != new

    def test_complete(self):
        for operator in ["order", "join"]:
            for strand in [1, -1]:
                self.check_ordering(strand=strand, operator=operator)


class TestConversion(unittest.TestCase):
    def test_convert_compound(self):
        parts = [
            BioSimple(1, 6, strand=1),
            BioSimple(10, 12, strand=1),
        ]
        bio = BioCompound(parts)
        location = CompoundLocation.from_biopython(bio)
        assert isinstance(location, CompoundLocation)
        assert not isinstance(bio, CompoundLocation)
        assert bio == location

    def test_convert_simple(self):
        for strand in [1, -1]:
            bio = BioSimple(1, 6, strand=strand)
            location = FeatureLocation.from_biopython(bio)
            assert isinstance(location, FeatureLocation)
            assert not isinstance(bio, FeatureLocation)
            assert bio == location


def get_distance_between_locations(first, second, **kwargs):
    static = _get_distance_between_locations(first, second, **kwargs)
    dynamic = first.get_distance_to(second, **kwargs)
    assert static == dynamic
    return static


class TestDistance(unittest.TestCase):
    def setUp(self):
        self.low = FeatureLocation(10, 20, 1)
        self.high = FeatureLocation(70, 80, 1)

    def test_no_wrapping(self):
        distance = get_distance_between_locations(self.low, self.high)
        assert distance == 50 == get_distance_between_locations(self.high, self.low)

    def test_wrapping(self):
        distance = get_distance_between_locations(self.low, self.high, wrap_point=100)
        assert distance == 30 == get_distance_between_locations(self.high, self.low, wrap_point=100)

    def test_closer_than_wrap(self):
        low = FeatureLocation(103, 3085, 1)
        high = FeatureLocation(3095, 4898, 1)
        wrap = 45016
        low_high = get_distance_between_locations(low, high, wrap_point=wrap)
        high_low = get_distance_between_locations(high, low, wrap_point=wrap)
        assert low_high == high_low
        assert low_high == 10

    def test_overlapping(self):
        self.high = FeatureLocation(15, 30, 1)
        assert locations_overlap(self.high, self.low)
        distance = get_distance_between_locations(self.low, self.high)
        assert distance == 0 == get_distance_between_locations(self.high, self.low)
