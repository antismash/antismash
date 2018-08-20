# License: GNU Affero General Public License v3 or later
# A copy of GNU AGPL v3 should have been included in this software package in LICENSE.txt.

# for test files, silence irrelevant and noisy pylint warnings
# pylint: disable=no-self-use,protected-access,missing-docstring

import unittest

from antismash.common.secmet import FeatureLocation
from antismash.common.secmet.features import SubRegion


class TestSubRegion(unittest.TestCase):
    def test_probability(self):
        loc = FeatureLocation(0, 10)
        assert SubRegion(loc, tool="test").probability is None
        assert SubRegion(loc, tool="test", probability=9.5).probability == 9.5
        assert SubRegion(loc, tool="test", probability=.5).probability == .5

    def test_anchor(self):
        loc = FeatureLocation(0, 10)
        assert SubRegion(loc, tool="test").anchor == ""
        assert SubRegion(loc, tool="test", anchor="anch").anchor == "anch"

    def test_orphaned(self):
        sub = SubRegion(FeatureLocation(0, 10), tool="test")
        assert not sub.parent_record
        with self.assertRaisesRegex(ValueError, "not in a record"):
            sub.get_subregion_number()
