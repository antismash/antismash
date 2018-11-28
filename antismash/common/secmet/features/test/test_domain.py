# License: GNU Affero General Public License v3 or later
# A copy of GNU AGPL v3 should have been included in this software package in LICENSE.txt.

# for test files, silence irrelevant and noisy pylint warnings
# pylint: disable=no-self-use,protected-access,missing-docstring

import unittest

from antismash.common.secmet import FeatureLocation
from antismash.common.secmet.features.domain import Domain


class TestDomain(unittest.TestCase):
    def test_construction(self):
        loc = FeatureLocation(1, 15, 1)
        domain = Domain(loc, "test_type", tool="test")
        assert domain.type == "test_type"
        assert domain.location == loc
        assert domain.created_by_antismash
        assert domain.tool == "test"
        assert domain.domain is None

    def test_translation(self):
        domain = Domain(FeatureLocation(1, 15, 1), "test_type", tool="test")
        with self.assertRaisesRegex(ValueError, "has no translation"):
            assert domain.translation is None
        domain.translation = "AAA"
        assert domain.translation == "AAA"

        with self.assertRaisesRegex(ValueError, "stop codons"):
            domain.translation = "A*A"

        for value in [7, None, Domain]:
            with self.assertRaises(AssertionError):
                domain.translation = value

        with self.assertRaisesRegex(ValueError, "empty"):
            domain.translation = ""
