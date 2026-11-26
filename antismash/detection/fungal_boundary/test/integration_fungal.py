# License: GNU Affero General Public License v3 or later
# A copy of GNU AGPL v3 should have been included in this software package in LICENSE.txt.

"""Integration test for fungal boundary detection on a saved Record."""

from typing import cast
import unittest
from unittest.mock import patch

from antismash.common.secmet import Record
from antismash.config import build_config, destroy_config
from antismash.custom_typing import AntismashModule
from antismash.detection import fungal_boundary
from antismash.detection.fungal_boundary.config import config


RECORD_PATH = (
    "antismash/detection/fungal_boundary/test/data/"
    "ANIA_07825_15genes_only_promoters.gbk"
)


def _run_saved_record() -> fungal_boundary.BoundaryDetectionResults:
    """Run the detection module on a subset of Aspergillus nidulans BN001304.1 contig."""
    record = Record.from_genbank(RECORD_PATH, taxon="fungi")[0]

    options = build_config(["--taxon", "fungi", "--boundaries-fungal"],
                           isolated=True, modules=[cast(AntismashModule, fungal_boundary)])
    try:
        return fungal_boundary.run_on_record(record, None, options)
    finally:
        destroy_config()


class TestFungalBoundaryIntegration(unittest.TestCase):
    """Check boundary predictions on a record"""

    def test_saved_record(self) -> None:
        """Check the reported anchors, boundaries and motifs."""
        mock_bg: dict[int, float] = {}
        for motif_id in range(config.n_motifs):
            mock_bg[motif_id] = 2e-3
        # The target motif (TCG-N5-CGA) gets a lower background frequency.
        mock_bg[869] = 5e-4

        # If running with this contig https://www.ncbi.nlm.nih.gov/nuccore/BN001304
        # which has the full aflatoxin BGC with ANIA_07825 as the core  
        # Then the predicted boundry should span from ANIA_11012 to ANIA_07825
        with patch.object(fungal_boundary.prediction, "get_motif_id_to_bg_freq", return_value=mock_bg):
            with patch.multiple(config, MIN_GENES_PER_CONTIG=15, MAX_GENES_IN_BGC=7):
                results = _run_saved_record()
        self.assertEqual(results.record_id, "BN001304.1")
        predictions = [
            (prediction.anchor_id, prediction.start, prediction.end, prediction.motif_string_dna,
             len(prediction.motif_counts))
            for prediction in results.predictions
        ]
        self.assertCountEqual(predictions, [
            ("ANIA_07825", 4870, 16337, "TCG-N5-CGA", 5),
        ])
