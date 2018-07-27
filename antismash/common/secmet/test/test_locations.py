# License: GNU Affero General Public License v3 or later
# A copy of GNU AGPL v3 should have been included in this software package in LICENSE.txt.

# for test files, silence irrelevant and noisy pylint warnings
# pylint: disable=no-self-use,protected-access,missing-docstring

import unittest

from antismash.common.secmet.locations import (
    convert_protein_position_to_dna,
    location_bridges_origin as is_bridged,
    split_origin_bridging_location as splitter,
    FeatureLocation,
    CompoundLocation
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

    def test_position_conversion_nonzero_start_compound(self):
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
