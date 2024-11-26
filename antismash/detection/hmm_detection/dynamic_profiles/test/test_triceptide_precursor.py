# License: GNU Affero General Public License v3 or later
# A copy of GNU AGPL v3 should have been included in this software package in LICENSE.txt.

# for test files, silence irrelevant and noisy pylint warnings
# pylint: disable=use-implicit-booleaness-not-comparison,protected-access,missing-docstring,consider-using-with

import unittest

from antismash.common.test.helpers import DummyRecord, DummyCDS

from antismash.detection.hmm_detection.dynamic_profiles import triceptide_precursor


def _create_cds(name: str, translation: str, offset: int = 0) -> tuple[DummyCDS, int]:
    end_coord = offset + len(translation) * 3
    return DummyCDS(start=offset, end=end_coord, strand=1, locus_tag=name, translation=translation), end_coord + 1


class TestTriceptidePrecursor(unittest.TestCase):
    def setUp(self) -> None:
        data = [
            ("YxD", "MTHTYLPGQPSEAQLPDFSALPLGELASCTEHPVLSTLLPGVCERLEDSGGAVAFYDDAPPVA"),
            ("HAA", "MMPLSPASAEPTASGTAVLDRVAARVRQRLETEQAATNRVGDGTHAASLIWPWPL"),
            ("YRD", "MLNKNHSVDQKMPINHPLLSRLIKEVENEQKTAMMFKYDRTHNRHNRGQ"),
            ("random_other_thing", "MAGICHATMAGICCATMAGICRAT"),
            ("too_short", "MAGIC"),
        ]

        features = []
        offset = 0
        for name, seq in data:
            feature, offset = _create_cds(name, seq, offset)
            features.append(feature)

        self.record = DummyRecord(features=features)

    def test_hits(self) -> None:
        results = triceptide_precursor.find_hits(self.record, {})
        assert set(results.keys()) == {"YxD", "HAA", "YRD"}
