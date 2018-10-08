# License: GNU Affero General Public License v3 or later
# A copy of GNU AGPL v3 should have been included in this software package in LICENSE.txt.

# for test files, silence irrelevant and noisy pylint warnings
# pylint: disable=no-self-use,protected-access,missing-docstring

import unittest

from antismash.common.test import helpers

from antismash.common.secmet.features.feature import (
    CompoundLocation,
    Feature,
    FeatureLocation,
    SeqFeature,
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
        # check that features without types are caught
        with self.assertRaises(AssertionError):
            sec = Feature.from_biopython(bio)
        bio.type = "test"
        sec = Feature.from_biopython(bio)
        assert sec.get_qualifier("foo") == tuple(["bar"])
        assert sec.get_qualifier("bar") is None

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

