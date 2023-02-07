# License: GNU Affero General Public License v3 or later
# A copy of GNU AGPL v3 should have been included in this software package in LICENSE.txt.

# for test files, silence irrelevant and noisy pylint warnings
# pylint: disable=use-implicit-booleaness-not-comparison,protected-access,missing-docstring

import unittest

from antismash.detection.cassis.islands import Island, get_islands
from antismash.detection.cassis.motifs import Motif
from antismash.detection.cassis.promoters import Promoter


class TestPairings(unittest.TestCase):
    def test_get_islands(self):
        motifs = [Motif(0, 3, hits={"gene1": 1, "gene2": 2}),
                  Motif(0, 4, hits={"gene2": 3, "gene4": 2, "gene5": 1})]
        # gene2 will be the anchor promoter
        anchor_promoter = 1
        promoters = []
        for i in range(1, 7):
            promoters.append(Promoter(f"gene{i}", i * 10, i * 10 + 4))
        # resulting in 2 different islands (this example)
        # promoter (pos): 1 2 3 4 5 6
        # binding sites:  1 2 0 0 0 0
        # island:         |-|
        first_island = Island(promoters[0], promoters[1], motifs[0])
        # promoter (pos): 1 2 3 4 5 6
        # binding sites:  0 3 0 2 1 0
        # island:           |---|
        second_island = Island(promoters[1], promoters[4], motifs[1])
        expected_islands = [first_island, second_island]
        assert get_islands(anchor_promoter, motifs, promoters) == expected_islands
