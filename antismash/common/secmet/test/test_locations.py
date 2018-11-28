# License: GNU Affero General Public License v3 or later
# A copy of GNU AGPL v3 should have been included in this software package in LICENSE.txt.

# for test files, silence irrelevant and noisy pylint warnings
# pylint: disable=no-self-use,protected-access,missing-docstring

import unittest

from antismash.common.secmet.locations import (
    convert_protein_position_to_dna,
    build_location_from_others,
    location_bridges_origin as is_bridged,
    split_origin_bridging_location as splitter,
    location_contains_other,
    location_from_string,
    locations_overlap,
    combine_locations,
    FeatureLocation,
    CompoundLocation,
    AfterPosition,
    BeforePosition,
    ExactPosition,
    UnknownPosition,
)


class TestProteinPositionConversion(unittest.TestCase):
    def setUp(self):
        self.func = convert_protein_position_to_dna

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

    def test_bad_strand(self):
        pairs = [(9, 12), (0, 3)]
        assert is_bridged(build_compound(pairs, 1))
        assert not is_bridged(build_compound(pairs, None))

    def test_not_bridged(self):
        assert not is_bridged(build_compound([(1, 6), (5, 10)], 1))
        assert not is_bridged(build_compound([(5, 10), (1, 6)], -1))


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

    def test_not_bridging_forward(self):
        loc = build_compound([(0, 3), (9, 12)], 1)
        with self.assertRaisesRegex(ValueError, "Location does not bridge origin"):
            print(splitter(loc))

    def test_not_bridging_reverse(self):
        loc = build_compound([(9, 12), (0, 3)], -1)
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


class TestCombiner(unittest.TestCase):
    def make(self, start, end):
        return FeatureLocation(start, end)

    def test_individual(self):
        loc = combine_locations(self.make(3, 7), self.make(5, 9))
        assert loc.start == 3 and loc.end == 9
        loc = combine_locations(self.make(3, 5), self.make(7, 9))
        assert loc.start == 3 and loc.end == 9
        loc = combine_locations(self.make(7, 9), self.make(3, 5))
        assert loc.start == 3 and loc.end == 9

        # it's silly, but since it theoretically is useful for CompoundLocation condensing
        loc = combine_locations(self.make(0, 5))
        assert loc.start == 0 and loc.end == 5
        loc = combine_locations(CompoundLocation([self.make(0, 3), self.make(6, 9)]))
        assert loc.start == 0 and loc.end == 9 and len(loc.parts) == 1

    def test_list(self):
        loc = combine_locations([self.make(i, i+1) for i in range(10, 20)])
        assert loc.start == 10 and loc.end == 20

    def test_generator(self):
        loc = combine_locations(i for i in [self.make(i, i+1) for i in range(10, 20)])
        assert loc.start == 10 and loc.end == 20

    def test_invalid(self):
        with self.assertRaisesRegex(TypeError, "object is not iterable"):
            combine_locations(0)
        with self.assertRaisesRegex(AttributeError, "has no attribute 'start'"):
            combine_locations(0, 1)


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
