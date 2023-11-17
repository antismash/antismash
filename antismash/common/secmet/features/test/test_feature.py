# License: GNU Affero General Public License v3 or later
# A copy of GNU AGPL v3 should have been included in this software package in LICENSE.txt.

# for test files, silence irrelevant and noisy pylint warnings
# pylint: disable=use-implicit-booleaness-not-comparison,protected-access,missing-docstring

import unittest
from unittest.mock import patch

from antismash.common.test import helpers

from antismash.common.secmet.features.feature import (
    CompoundLocation,
    Feature,
    FeatureLocation,
    pop_locus_qualifier as pop_locus,
    SeqFeature,
)
from antismash.common.secmet import features
from antismash.common.secmet.locations import (
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
        with self.assertRaisesRegex(ValueError, "invalid length"):
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
        for feature_type in ["protocluster", "cds_motif", "test"]:
            for start, end in [(1, 5), (3, 8), (10, 15)]:
                feature = Feature(FeatureLocation(start, end, strand=1),
                                  feature_type=feature_type)
                assert str(feature) == f"{feature_type}([{start}:{end}](+))"
                feature = Feature(FeatureLocation(start, end, strand=-1),
                                  feature_type=feature_type)
                assert str(feature) == f"{feature_type}([{start}:{end}](-))"

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
        for codon_start in ["-1", "4", "NA"]:
            seqf.qualifiers["codon_start"] = [codon_start]
            with self.assertRaisesRegex(ValueError, "invalid codon_start"):
                Feature.from_biopython(seqf)

    def test_feature_type(self):
        loc = FeatureLocation(1, 4, strand=1)
        with self.assertRaisesRegex(ValueError, "invalid length"):
            Feature(loc, feature_type="")
        assert Feature(loc, feature_type="A")
        assert Feature(loc, feature_type="A"*15)
        with self.assertRaisesRegex(ValueError, "invalid length"):
            Feature(loc, feature_type="A"*16)

    def test_ordering(self):
        feature = Feature(FeatureLocation(10, 40), feature_type="test")
        before = Feature(FeatureLocation(5, 30), feature_type="test")
        after = Feature(FeatureLocation(20, 40), feature_type="test")
        longer = Feature(FeatureLocation(10, 70), feature_type="test")

        def check(first, second):
            assert first < second
            assert first < second.location
            assert sorted([second, first]) == [first, second]

        check(before, feature)
        check(feature, after)
        check(feature, longer)

        assert sorted([feature, before, after, longer]) == [before, feature, longer, after]

    def test_negative_locations(self):
        loc = FeatureLocation(-2, 500, strand=1)
        with self.assertRaisesRegex(ValueError, "negative coordinate"):
            Feature(loc, feature_type="")


class TestSimpleCoordinates(unittest.TestCase):
    def test_simple_location(self):
        for strand in [1, None, -1]:
            feature = Feature(FeatureLocation(5, 10, strand), "test")
            assert feature.start == 5
            assert feature.end == 10

    def test_compound_location(self):
        for strand in [1, None, -1]:
            parts = [FeatureLocation(5, 10, strand), FeatureLocation(15, 20, strand)]
            if strand is not None:
                parts = parts[::strand]
            feature = Feature(CompoundLocation(parts), "test")
            assert feature.start == 5
            assert feature.end == 20

    def test_cross_origin_location(self):
        for strand in [1, None, -1]:
            parts = [FeatureLocation(80, 100, strand), FeatureLocation(5, 20, strand)]
            if strand is not None:
                parts = parts[::strand]
            feature = Feature(CompoundLocation(parts), "test")
            assert feature.start == 80
            assert feature.end == 20


class TestSubLocation(unittest.TestCase):
    def setUp(self):
        self.feature = Feature(FeatureLocation(10, 40, 1), feature_type="test")
        self.get_sub = self.feature.get_sub_location_from_protein_coordinates

    def test_invalid(self):
        for bad_start, bad_end in [(-1, 1), (1, -1), (1, 11)]:
            with self.assertRaisesRegex(ValueError, "must be contained by the feature"):
                self.get_sub(bad_start, bad_end)
        for bad_start, bad_end in [("test", 5), (5, "test"), (None, 5)]:
            with self.assertRaisesRegex(TypeError, "(unorderable types|not supported)"):
                self.get_sub(bad_start, bad_end)
        with self.assertRaisesRegex(ValueError, "must be less than the end"):
            self.get_sub(5, 1)

    @patch.object(features.feature, "convert_protein_position_to_dna", return_value=(9, 15))
    def test_start_outside(self, _mocked):
        with self.assertRaisesRegex(ValueError, "Protein coordinate start .* is outside feature"):
            self.get_sub(1, 5)

    @patch.object(features.feature, "convert_protein_position_to_dna", return_value=(15, 41))
    def test_end_outside(self, _mocked):
        with self.assertRaisesRegex(ValueError, "Protein coordinate end .* is outside feature"):
            self.get_sub(1, 5)

    @patch.object(features.feature, "convert_protein_position_to_dna", return_value=(10, 3))
    def test_inverted(self, _mocked):
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


class TestLocusPop(unittest.TestCase):
    def test_base(self):
        for value in ["locus", "something"]:
            assert pop_locus({"locus_tag": [value]}) == value

    def test_default(self):
        for default in ["test", None]:
            assert pop_locus({}, default=default) == default

    def test_required(self):
        with self.assertRaises(KeyError):
            pop_locus({}, allow_missing=False)
        with self.assertRaises(KeyError):
            pop_locus({}, allow_missing=False, default="some value")

    def test_biopython_spaces(self):
        locus = "longnamewith space"
        assert pop_locus({"locus_tag": [locus]}) == locus.replace(" ", "")
