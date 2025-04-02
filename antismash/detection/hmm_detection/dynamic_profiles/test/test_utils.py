# License: GNU Affero General Public License v3 or later
# A copy of GNU AGPL v3 should have been included in this software package in LICENSE.txt.

# for test files, silence irrelevant and noisy pylint warnings
# pylint: disable=use-implicit-booleaness-not-comparison,protected-access,missing-docstring,consider-using-with

import re
import unittest

from antismash.common.test.helpers import DummyRecord, DummyCDS

from antismash.detection.hmm_detection.dynamic_profiles import _utils as utils


class TestFindMotifAroundAnchor(unittest.TestCase):
    def setUp(self) -> None:
        #      0    5    10   15   20   25   30   35   40   45   50   65   60   65   70   75   80   85   90   95   100  105  110  115
        #            M  A  G  I  C  H  A  T  *      M  A  I  N  G  E  N  E  *         M  I  C  H  *        M  A  G  I  C  H  A  T  *
        seq = "AAAAAATGGCGGGGATTTGTCACGCGACGTGAAAAAATGGCGATTAACGGGGAGAACGAGTGAAAAAAAAATGATTTGTCACTGAAAAAAAATGGCGGGGATTTGTCACGCGACGTGA"
        self.anchor = DummyCDS(36, 62, locus_tag="MAIN", translation="MAINGENE")
        self.record = DummyRecord(features=[self.anchor], seq=seq, record_id="Dummy")
        self.motif = re.compile("ICH")

    def run_func(self, max_dist, min_len, max_len=0, early_abort=False):
        return utils.find_motif_around_anchor(self.record, self.anchor, self.motif,
                                              max_dist=max_dist, min_len=min_len,
                                              max_len=max_len, early_abort=early_abort)

    def test_simple(self):
        found = self.run_func(max_dist=100, min_len=3)
        assert len(found) == 3, f"{found} has wrong length"
        assert found[0].location.start == 5
        assert found[1].location.start == 70
        assert found[2].location.start == 91

    def test_max_dist(self):
        found = self.run_func(max_dist=50, min_len=3)
        assert len(found) == 2, f"{found} has wrong length"
        assert found[0].location.start == 5
        assert found[1].location.start == 70

    def test_min_len(self):
        found = self.run_func(max_dist=100, min_len=5)
        assert len(found) == 2, f"{found} has wrong length"
        assert found[0].location.start == 5
        assert found[1].location.start == 91

    def test_max_len(self):
        found = self.run_func(max_dist=100, min_len=3, max_len=5)
        assert len(found) == 1, f"{found} has wrong length"
        assert found[0].location.start == 70

    def test_early_abort(self):
        found = self.run_func(max_dist=100, min_len=3, early_abort=True)
        assert len(found) == 1, f"{found} has wrong length"
        assert found[0].location.start == 5

    def test_raises_max_len(self):
        with self.assertRaisesRegex(ValueError, "max_len needs to be >= 0"):
            self.run_func(max_dist=100, min_len=3, max_len=-1)


class TestFindMotifAroundAnchorAnnotated(TestFindMotifAroundAnchor):
    def setUp(self) -> None:
        super().setUp()
        self.record.add_feature(DummyCDS(5, 31, locus_tag="FIRST", translation="MAGICHAT"))
        self.record.add_feature(DummyCDS(70, 84, locus_tag="SMALL", translation="MICH"))
        self.record.add_feature(DummyCDS(91, 117, locus_tag="SECON", translation="MAGICHAT"))

