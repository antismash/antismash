# License: GNU Affero General Public License v3 or later
# A copy of GNU AGPL v3 should have been included in this software package in LICENSE.txt.

# for test files, silence irrelevant and noisy pylint warnings
# pylint: disable=use-implicit-booleaness-not-comparison,protected-access,missing-docstring,consider-using-with

import unittest
from unittest.mock import MagicMock

from antismash.common.test.helpers import DummyRecord, DummyCDS

from antismash.detection.hmm_detection.dynamic_profiles import triceptide_rsam


def _create_cds(name: str, translation: str, offset: int = 0) -> tuple[DummyCDS, int]:
    end_coord = offset + len(translation) * 3
    return DummyCDS(start=offset, end=end_coord, strand=1, locus_tag=name, translation=translation), end_coord + 1


class TestTriceptideRsam(unittest.TestCase):
    def setUp(self) -> None:

        data = [
            ("rSAM", "MAGIC"*20),
            ("HAA", "MMPLSPASAEPTASGTAVLDRVAARVRQRLETEQAATNRVGDGTHAASLIWPWPL"),
        ]

        features = []
        offset = 0
        for name, seq in data:
            feature, offset = _create_cds(name, seq, offset)
            features.append(feature)

        self.record = DummyRecord(features=features)

    def test_hits_single_profile(self) -> None:
        results = triceptide_rsam.find_hits(self.record, {"rSAM": [MagicMock(query_id="PTHR43273")]})
        assert set(results.keys()) == {"rSAM"}

    def test_hits_two_domains(self) -> None:
        results = triceptide_rsam.find_hits(self.record, {"rSAM": [MagicMock(query_id="PF04055"), MagicMock(query_id="SPASM")]})
        assert set(results.keys()) == {"rSAM"}

    def test_hits_one_missing(self) -> None:
        results = triceptide_rsam.find_hits(self.record, {"rSAM": [MagicMock(query_id="PF04055")]})
        assert not set(results.keys())

    def test_hits_no_hits(self) -> None:
        results = triceptide_rsam.find_hits(self.record, {})
        assert not set(results.keys())
