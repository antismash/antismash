# License: GNU Affero General Public License v3 or later
# A copy of GNU AGPL v3 should have been included in this software package in LICENSE.txt.

# for test files, silence irrelevant and noisy pylint warnings
# pylint: disable=use-implicit-booleaness-not-comparison,protected-access,missing-docstring

import unittest

from antismash.common.secmet.locations import CompoundLocation, FeatureLocation
from antismash.support.genefinding.run_prodigal import (
    _build_location_from_prodigal,
    ProdigalGene,
)


def build_location(start, end, strand, length):
    gene = ProdigalGene("dummy", start, end, strand)
    return _build_location_from_prodigal(gene, length)


class TestProdigal(unittest.TestCase):
    def test_build_compound(self):
        length = 20
        for strand in [1, -1]:
            loc = build_location(15, 5, strand, length)  # 1-indexed start
            print(strand, loc)
            assert loc == CompoundLocation([
                FeatureLocation(15, length, strand),
                FeatureLocation(0, 5, strand),
            ][::strand])

    def test_simple(self):
        length = 20
        for strand in [1, -1]:
            loc = build_location(5, 15, strand, length)  # 1-indexed start
            assert loc == FeatureLocation(5, 15, strand)

    def test_bad_lengths(self):
        for length in [-1, 5, 14]:  # anything less than the start when start < end
            with self.assertRaisesRegex(ValueError, "exceeds max length"):
                build_location(15, 5, 1, length)

    def test_edge(self):
        loc = build_location(5, 15, 1, length=15)  # 1-indexed start
        assert loc == FeatureLocation(5, 15, 1)
