# License: GNU Affero General Public License v3 or later
# A copy of GNU AGPL v3 should have been included in this software package in LICENSE.txt.

# for test files, silence irrelevant and noisy pylint warnings
# pylint: disable=no-self-use,protected-access,missing-docstring

import unittest

from minimock import mock, restore

from antismash.common.test import helpers

from antismash.common.secmet.features.feature import (
    CompoundLocation,
    Feature,
    FeatureLocation,
    SeqFeature,
    _adjust_location_by_offset as adjust,
)
from antismash.common.secmet import features  # mocked, pylint: disable=unused-import
from antismash.common.secmet.locations import (
    ExactPosition,
    BeforePosition,
    AfterPosition,
)


class TestFeature(unittest.TestCase):
    def test_overlaps_with(self):
        # no overlap
        feature = helpers.DummyFeature(5, 10)
        assert isinstance(feature, Feature)  # just to be sure it works the way we want
        other = helpers.DummyFeature(100, 110)
        assert not feature.overlaps_with(other) and not other.overlaps_with(feature)
        # completely within
        other = helpers.DummyFeature(0, 20)
        assert feature.overlaps_with(other) and other.overlaps_with(feature)
        # partially within
        other = helpers.DummyFeature(0, 8)
        assert feature.overlaps_with(other) and other.overlaps_with(feature)
        # borders touching
        other = helpers.DummyFeature(0, 5)
        assert not (feature.overlaps_with(other) or other.overlaps_with(feature))

    def test_is_contained_by(self):
        # same location is considered to be contained
        feature = helpers.DummyFeature(5, 10)
        assert feature.is_contained_by(feature)
        for strand in (-1, 1):
            # no overlap
            other = helpers.DummyFeature(15, 25, strand)
            assert not feature.is_contained_by(other)
            assert not other.is_contained_by(feature)
            # b is contained
            other = helpers.DummyFeature(6, 9, strand)
            assert not feature.is_contained_by(other)
            assert other.is_contained_by(feature)
            # only partial overlap
            other = helpers.DummyFeature(6, 19, strand)
            assert not feature.is_contained_by(other)
            assert not other.is_contained_by(feature)
            other = helpers.DummyFeature(1, 7, strand)
            assert not feature.is_contained_by(other)
            assert not other.is_contained_by(feature)
            # edge cases
            other = helpers.DummyFeature(5, 7, strand)
            assert not feature.is_contained_by(other)
            assert other.is_contained_by(feature)
            other = helpers.DummyFeature(7, 10, strand)
            assert not feature.is_contained_by(other)
            assert other.is_contained_by(feature)

    def test_biopython_conversion(self):
        bio = SeqFeature(FeatureLocation(1, 5))
        bio.qualifiers["foo"] = ["bar"]
        bio.qualifiers["key_only"] = None
        # check that features without types are caught
        with self.assertRaises(AssertionError):
            sec = Feature.from_biopython(bio)
        bio.type = "test"
        sec = Feature.from_biopython(bio)
        assert sec.get_qualifier("foo") == tuple(["bar"])
        assert sec.get_qualifier("bar") is None
        assert sec.get_qualifier("key_only") is True

    def test_created_by_antismash_conversion(self):
        for created in [True, False]:
            old = Feature(FeatureLocation(1, 5), feature_type="testtype", created_by_antismash=created)
            assert old.created_by_antismash == created
            assert old.type == "testtype"

            new = Feature.from_biopython(old.to_biopython()[0])
            assert new.created_by_antismash == old.created_by_antismash == created
            assert new.type == old.type == "testtype"

    def test_string_conversion(self):
        for feature_type in ["cluster", "cds_motif", "test"]:
            for start, end in [(1, 5), (3, 8), (10, 15)]:
                feature = Feature(FeatureLocation(start, end, strand=1),
                                  feature_type=feature_type)
                assert str(feature) == "%s([%d:%d](+))" % (feature_type, start, end)
                feature = Feature(FeatureLocation(start, end, strand=-1),
                                  feature_type=feature_type)
                assert str(feature) == "%s([%d:%d](-))" % (feature_type, start, end)

    def test_bridging_fails(self):
        parts = [FeatureLocation(9, 12, strand=1), FeatureLocation(0, 3, strand=1)]
        with self.assertRaisesRegex(ValueError, "bridge the record origin"):
            Feature(CompoundLocation(parts, operator="join"), feature_type="test")
        Feature(CompoundLocation(parts[::-1], operator="join"), feature_type="test")

    def test_conversion_with_codon_start_forward(self):
        seqf = SeqFeature(FeatureLocation(BeforePosition(5), 12), 1)
        seqf.type = "test"
        for codon_start in "123":
            seqf.qualifiers["codon_start"] = [codon_start]
            feature = Feature.from_biopython(seqf)
            assert feature._original_codon_start == int(codon_start) - 1
            assert feature.location.start == BeforePosition(5 + feature._original_codon_start)
            new = feature.to_biopython()[0]
            assert new.qualifiers["codon_start"] == [codon_start]
            assert new.location.start == seqf.location.start
            assert new.location.end == seqf.location.end

    def test_conversion_with_codon_start_reverse(self):
        seqf = SeqFeature(FeatureLocation(5, AfterPosition(12), -1))
        seqf.type = "test"
        for codon_start in "123":
            seqf.qualifiers["codon_start"] = [codon_start]
            feature = Feature.from_biopython(seqf)
            assert feature._original_codon_start == int(codon_start) - 1
            assert feature.location.end == AfterPosition(12 - feature._original_codon_start)
            new = feature.to_biopython()[0]
            assert new.qualifiers["codon_start"] == [codon_start]
            assert new.location.start == seqf.location.start
            assert new.location.end == seqf.location.end

    def test_invalid_codon_start(self):
        seqf = SeqFeature(FeatureLocation(5, AfterPosition(12), -1))
        seqf.type = "test"
        for codon_start in ["-1", "4"]:
            seqf.qualifiers["codon_start"] = [codon_start]
            with self.assertRaisesRegex(ValueError, "invalid codon_start"):
                Feature.from_biopython(seqf)


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


class TestSubLocation(unittest.TestCase):
    def setUp(self):
        self.feature = Feature(FeatureLocation(10, 40, 1), feature_type="test")
        self.get_sub = self.feature.get_sub_location_from_protein_coordinates

    def tearDown(self):
        restore()

    def test_invalid(self):
        for bad_start, bad_end in [(-1, 1), (1, -1), (1, 11)]:
            with self.assertRaisesRegex(ValueError, "must be contained by the feature"):
                self.get_sub(bad_start, bad_end)
        for bad_start, bad_end in [("test", 5), (5, "test"), (None, 5)]:
            with self.assertRaisesRegex(TypeError, "(unorderable types|not supported)"):
                self.get_sub(bad_start, bad_end)
        with self.assertRaisesRegex(ValueError, "must be less than the end"):
            self.get_sub(5, 1)
        mock("features.feature.convert_protein_position_to_dna", returns=(9, 15))
        with self.assertRaisesRegex(ValueError, "Protein coordinate start .* is outside feature"):
            self.get_sub(1, 5)
        mock("features.feature.convert_protein_position_to_dna", returns=(15, 41))
        with self.assertRaisesRegex(ValueError, "Protein coordinate end .* is outside feature"):
            self.get_sub(1, 5)
        mock("features.feature.convert_protein_position_to_dna", returns=(10, 3))
        with self.assertRaisesRegex(ValueError, "Invalid protein coordinate conversion"):
            self.get_sub(1, 5)

    def test_simple_forward(self):
        assert self.get_sub(0, 1) == FeatureLocation(10, 13, 1)
        assert self.get_sub(2, 4) == FeatureLocation(16, 22, 1)
        assert self.get_sub(9, 10) == FeatureLocation(37, 40, 1)

    def test_simple_reverse(self):
        self.feature.location = FeatureLocation(10, 40, -1)
        assert self.get_sub(0, 1) == FeatureLocation(37, 40, -1)
        assert self.get_sub(2, 4) == FeatureLocation(28, 34, -1)
        assert self.get_sub(9, 10) == FeatureLocation(10, 13, -1)

    def test_compound_reverse(self):
        self.feature.location = CompoundLocation([FeatureLocation(21, 27, -1),
                                                  FeatureLocation(12, 15, -1),
                                                  FeatureLocation(0, 6, -1)])
        assert self.get_sub(2, 3) == FeatureLocation(12, 15, -1)

    def test_compound_overlap_forward(self):
        self.feature.location = CompoundLocation([FeatureLocation(10, 16, 1),
                                                  FeatureLocation(15, 24, 1)])
        assert self.get_sub(0, 1) == FeatureLocation(10, 13, 1)
        assert self.get_sub(1, 3) == CompoundLocation([FeatureLocation(13, 16, 1),
                                                       FeatureLocation(15, 18, 1)])
        assert self.get_sub(4, 5) == FeatureLocation(21, 24, 1)

    def test_compound_overlap_reverse(self):
        self.feature.location = CompoundLocation([FeatureLocation(15, 24, -1),
                                                  FeatureLocation(10, 16, -1)])
        assert self.get_sub(0, 1) == FeatureLocation(21, 24, -1)
        assert self.get_sub(2, 4) == CompoundLocation([FeatureLocation(15, 18, -1),
                                                       FeatureLocation(13, 16, -1)])
        assert self.get_sub(4, 5) == FeatureLocation(10, 13, -1)
