# License: GNU Affero General Public License v3 or later
# A copy of GNU AGPL v3 should have been included in this software package in
# LICENSE.txt.

"""Unit tests for fungal boundary prediction utils script"""

# pylint: disable=protected-access

import unittest

import numpy as np

from .. import utils
from ..config import config


class TestMotifCounts(unittest.TestCase):
    """Test selection of nonzero motif counts and count capping."""

    def test_capped_counts(self) -> None:
        """Check that nonzero motif IDs are retained and counts are capped."""
        counts = np.array([0, 2, 0, 3, 5], dtype=np.int64)

        result = utils.MotifCounts.from_counts(counts, capped_max_count=3)

        np.testing.assert_array_equal(result.motif_ids, [1, 3, 4])
        np.testing.assert_array_equal(result.counts, [2, 3, 3])


class TestGetMotifType(unittest.TestCase):
    """Check motif types at and just below ID cutoffs"""

    def test_default_boundaries(self) -> None:
        """Check the IDs on either side of both default boundaries."""
        direct, linear = config.motif_type_id_boundaries

        assert utils.get_motif_type(direct - 1) is utils.MotifType.INVERTED
        assert utils.get_motif_type(direct) is utils.MotifType.DIRECT
        assert utils.get_motif_type(linear - 1) is utils.MotifType.DIRECT
        assert utils.get_motif_type(linear) is utils.MotifType.LINEAR


class TestMotifText(unittest.TestCase):
    """Test conversion of motif IDs to DNA patterns."""

    def test_inverted_repeat(self) -> None:
        """Check that the second arm is the reverse complement."""
        motif_id = 6 * config.spacers + 5  # ACG with a five-base spacer
        assert utils.convert_motif_id_to_text(motif_id) == "ACG-N5-CGT"

    def test_direct_repeat(self) -> None:
        """Check that both arms are identical."""
        motif_id = config.direct_motif_id_start + 6 * config.spacers + 5
        assert utils.convert_motif_id_to_text(motif_id) == "ACG-N5-ACG"

    def test_linear_motif(self) -> None:
        """Check that decoding preserves the full six-base length."""
        assert utils.convert_motif_id_to_text(config.linear_motif_id_start) == "AAAAAA"


class TestKmerEncoding(unittest.TestCase):
    """Check that k-mer IDs, reverse-complement IDs, and validity flags follows logic"""

    def test_overlapping_kmers(self) -> None:
        """Encode AAA, AAC, and ACG in order, including the valid ID zero."""
        encoded, valid, complements = utils._encode_kmers("AAACG", 3)

        np.testing.assert_array_equal(encoded, [0, 1, 6])
        np.testing.assert_array_equal(valid, [True, True, True])
        np.testing.assert_array_equal(complements, [63, 47, 27])

    def test_ambiguous_base(self) -> None:
        """Flag windows containing N while keeping their positions."""
        encoded, valid, complements = utils._encode_kmers("AACNGTA", 3)

        np.testing.assert_array_equal(valid, [True, False, False, False, True])
        np.testing.assert_array_equal(encoded[valid], [1, 44])
        np.testing.assert_array_equal(complements[valid], [47, 49])

    def test_sequence_shorter_than_k(self) -> None:
        """Return three empty arrays when a full window cannot fit."""
        for sequence in ("", "AC"):
            with self.subTest(sequence=sequence):
                encoded, valid, complements = utils._encode_kmers(sequence, 3)
                assert encoded.shape == valid.shape == complements.shape == (0,)

    def test_sequence_length_equals_k(self) -> None:
        """Encode exactly one window when the sequence length equals k."""
        encoded, valid, complements = utils._encode_kmers("AC", 2)

        np.testing.assert_array_equal(encoded, [1])
        np.testing.assert_array_equal(valid, [True])
        np.testing.assert_array_equal(complements, [11])

    def test_reverse_complement_lookup(self) -> None:
        """Check the IDs and both directions of a reverse complement pair."""
        kmer_ids, complement_ids = utils._get_kmer_to_ids(3)

        assert kmer_ids["ACG"] == 6
        assert kmer_ids["CGT"] == 27
        assert complement_ids[6] == 27
        assert complement_ids[27] == 6


class TestRepeatMotifs(unittest.TestCase):
    """Test repeat motif counting, spacers, and N bases."""

    def test_inverted_repeat(self) -> None:
        """Count reverse-complement arms with -and without a spacer."""
        base_motif_id = 6 * config.spacers  # ACG with no spacer
        for spacer in ("", "CC"):
            with self.subTest(spacer=spacer):
                result = utils.count_repeat_motifs(f"ACG{spacer}CGT")
                np.testing.assert_array_equal(result.motif_ids, [base_motif_id + len(spacer)])
                np.testing.assert_array_equal(result.counts, [1])

    def test_direct_repeat(self) -> None:
        """Count identical arms with -and without a spacer."""
        arm_id = 6  # ACG
        base_motif_id = config.direct_motif_id_start + arm_id * config.spacers
        for spacer in ("", "CC"):
            with self.subTest(spacer=spacer):
                result = utils.count_repeat_motifs(f"ACG{spacer}ACG")
                np.testing.assert_array_equal(result.motif_ids, [base_motif_id + len(spacer)])
                np.testing.assert_array_equal(result.counts, [1])


class TestLinearMotifs(unittest.TestCase):
    """Test linear motif counting and reverse-complement grouping."""

    def test_reverse_complements_share_motif_id(self) -> None:
        """Count both orientations under the same motif ID."""
        result = utils.count_linear_motifs("AAAAAANTTTTTT", capped_max_count=None)

        np.testing.assert_array_equal(result.motif_ids, [config.linear_motif_id_start])
        np.testing.assert_array_equal(result.counts, [2])
