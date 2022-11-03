# License: GNU Affero General Public License v3 or later
# A copy of GNU AGPL v3 should have been included in this software package in LICENSE.txt.

# for test files, silence irrelevant and noisy pylint warnings
# pylint: disable=use-implicit-booleaness-not-comparison,protected-access,missing-docstring

import unittest

from antismash.detection.cassis.pairings import Pairing, PROMOTER_RANGE


class TestPairings(unittest.TestCase):
    def test_pairing_count(self):
        self.assertEqual(len(PROMOTER_RANGE), 250)

    def test_pairing_types(self):
        assert isinstance(PROMOTER_RANGE[0], Pairing)

    def test_pairing_string(self):
        assert PROMOTER_RANGE[0].pairing_string == "+00_-03"
