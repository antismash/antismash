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
    "ANIA_07825_300genes_only_promoters.gbk"
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
        with patch.object(config, "MIN_GENES_PER_CONTIG", 300):
            results = _run_saved_record()
        self.assertEqual(results.record_id, "BN001304.1")
        predictions = [
            (prediction.anchor_id, prediction.start, prediction.end, prediction.motif_string_dna)
            for prediction in results.predictions
        ]
        self.assertCountEqual(predictions, [
            ("ANIA_07825", 226599, 246313, "TCG-N5-CGA"),
        ])
