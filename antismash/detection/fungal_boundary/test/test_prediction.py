# License: GNU Affero General Public License v3 or later
# A copy of GNU AGPL v3 should have been included in this software package in
# LICENSE.txt.

"""tests for fungal boundary prediction."""

import unittest
from unittest.mock import Mock

import numpy as np

from .. import prediction


def gene_with_promoter_counts(motif_ids: list[int], counts: list[int]) -> prediction.GeneWithMotifs:
    """Create a test gene with the given promoter motif IDs and counts."""
    return Mock(
        spec=prediction.GeneWithMotifs,
        motif_ids=np.array(motif_ids),
        capped_motif_counts=np.array(counts),
    )


class TestBackgroundFrequencies(unittest.TestCase):
    """Test background frequencies for observed motifs."""

    def test_background_frequencies_by_motif_type(self) -> None:
        """Sum promoter counts across genes and normalise within each motif type."""
        genes = {
            "A": gene_with_promoter_counts([0, 10], counts=[3, 1]),
            "B": gene_with_promoter_counts([0, 1, 20], counts=[1, 2, 1]),
        }

        # to test the different motif type (inverted, etc) motif_type_id_boundaries is defined. 
        result = prediction.get_motif_id_to_bg_freq(genes, motif_type_id_boundaries=(10, 20))

        assert result == {0: 0.75, 1: 0.25, 10: 1.0, 20: 1.0}
