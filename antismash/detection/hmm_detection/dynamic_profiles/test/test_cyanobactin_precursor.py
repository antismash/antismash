import unittest

from antismash.common.hmm_rule_parser.structures import DynamicHit, DynamicProfile
from antismash.common.secmet import Record
from antismash.common.test.helpers import DummyRecord, DummyCDS

from antismash.detection.hmm_detection.dynamic_profiles import cyanobactin_precursor


def _create_cds(name: str, translation: str, offset: int = 0) -> tuple[DummyCDS, int]:
    end_coord = offset + len(translation) * 3
    return DummyCDS(start=offset, end=end_coord, strand=1, locus_tag=name, translation=translation), end_coord + 1


class TestCyanobactinPrecursor(unittest.TestCase):
    def setUp(self) -> None:
        data = [
            ("exact_match", "MTKKNIRPQQVAPVERETISTAKDQSGQVQAQTSQIWGSPVPFAGDDAE"),
            ("exact_match_with_offset",
             "MRITPMDKKNILPQQGKPVFRTTTGKLPSYLAELSEEALGGNGLEASHCATICAFDGAEASHCATICAFDGDEA"),
            ("one_mismatch", "MTAKNIRPQQVAPVERETISTAKDQSGQVQAQTSQIWGSPVPFAGDDAE"),
            ("two_mismatches", "MTAANIRPQQVAPVERETISTAKDQSGQVQAQTSQIWGSPVPFAGDDAE"),
            ("random_other_thing", "MAGICHATMAGICCATMAGICRAT"),
            ("too_short", "MAGIC"),
            ("off_by_one_check", "MCKTLHDTNEGKNLTPFSSGPVR"),
            ("pv_too_early", "MPVMAGICHATMAGICCATMAGICRAT"),
            ("kk_too_late", "MAGICHATMAGICCATMAGICRATKKN")
        ]

        features = []
        offset = 0
        for name, seq in data:
            feature, offset = _create_cds(name, seq, offset)
            features.append(feature)

        self.record = DummyRecord(features=features)

    def test_one_mismatch(self) -> None:
        results = cyanobactin_precursor.find_hits(self.record)
        assert set(results.keys()) == {"exact_match", "one_mismatch", "exact_match_with_offset"}

    def test_two_mismatches(self) -> None:
        results = cyanobactin_precursor.find_hits(self.record, 2)
        assert set(results.keys()) == {"exact_match", "one_mismatch", "two_mismatches", "exact_match_with_offset"}
