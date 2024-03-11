# License: GNU Affero General Public License v3 or later
# A copy of GNU AGPL v3 should have been included in this software package in LICENSE.txt.

# for test files, silence irrelevant and noisy pylint warnings
# pylint: disable=use-implicit-booleaness-not-comparison,protected-access,missing-docstring

import unittest

from antismash.common.secmet.locations import FeatureLocation
from antismash.support.genefinding.run_prodigal import (
    _build_location_from_prodigal,
    ProdigalGene,
)


def build_location(start, end, strand):
    gene = ProdigalGene("dummy", start, end, strand)
    return _build_location_from_prodigal(gene)


class TestProdigal(unittest.TestCase):
    def test_simple(self):
        for strand in [1, -1]:
            loc = build_location(5, 15, strand)  # 1-indexed start
            assert loc == FeatureLocation(5, 15, strand)
