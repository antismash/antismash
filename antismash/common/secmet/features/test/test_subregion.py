# License: GNU Affero General Public License v3 or later
# A copy of GNU AGPL v3 should have been included in this software package in LICENSE.txt.

# for test files, silence irrelevant and noisy pylint warnings
# pylint: disable=no-self-use,protected-access,missing-docstring

import unittest

from antismash.common.secmet import FeatureLocation
from antismash.common.secmet.features.subregion import SubRegion, SideloadedSubRegion


class TestSubRegion(unittest.TestCase):
    def test_anchor(self):
        loc = FeatureLocation(0, 10)
        assert SubRegion(loc, tool="test").label == ""
        assert SubRegion(loc, tool="test", label="anch").label == "anch"

    def test_orphaned(self):
        sub = SubRegion(FeatureLocation(0, 10), tool="test")
        assert not sub.parent_record
        with self.assertRaisesRegex(ValueError, "not in a record"):
            sub.get_subregion_number()


class TestSideloaded(unittest.TestCase):
    def test_conversion(self):
        loc = FeatureLocation(10, 20)
        extras = {"a": ["5", "c"], "b": ["something"]}
        source = SideloadedSubRegion(loc, "tool name", label="some label", extra_qualifiers=extras)

        bio_features = source.to_biopython()
        assert len(bio_features) == 1
        for key, val in extras.items():
            assert bio_features[0].qualifiers[key] == val

        for regenerator in [SideloadedSubRegion, SubRegion]:
            dest = regenerator.from_biopython(bio_features[0])
            assert isinstance(dest, SideloadedSubRegion)
            assert dest.extra_qualifiers == source.extra_qualifiers == extras
            assert dest.tool == source.tool
            assert dest.label == source.label
            assert dest.location == source.location

            for key, val in extras.items():
                assert not dest.get_qualifier(key)
