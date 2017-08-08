# License: GNU Affero General Public License v3 or later
# A copy of GNU AGPL v3 should have been included in this software package in LICENSE.txt.

import unittest

from Bio.SeqFeature import BeforePosition, AfterPosition, ExactPosition

from antismash.common.test.helpers import FakeRecord
from antismash.modules.genefinding.all_orfs import find_all_orfs, scan_orfs

class TestOrfCounts(unittest.TestCase):
    def run_both_dirs(self, expected, seq):
        assert expected == find_all_orfs(FakeRecord(seq=seq, real_seq=True))

        # replace for reverse direction testing
        for start, comp in [('ATG', 'TAC'), ('GTG', 'CAC'), ('TTG', 'AAC')]:
            seq = seq.replace(start, comp)
        for stop, comp in [('TAA', 'ATT'), ('TAG', 'ATC'), ('TGA', 'ACT')]:
            seq = seq.replace(stop, comp)
        assert expected == find_all_orfs(FakeRecord(seq=seq, real_seq=True))

    def test_empty_sequence(self):
        self.run_both_dirs(0, "")

    def test_no_hits(self):
        # no starts or stops
        self.run_both_dirs(0, "XXX")

        # nothing > 60
        self.run_both_dirs(0, "ATGXXXTGA")
        self.run_both_dirs(0, "ATG"+ "X"*54 +"TGA")

    def test_all_combos(self):
        for start in ('ATG', 'GTG', 'TTG'):
            for stop in ('TAA', 'TAG', 'TGA'):
                seq = "{}{}{}".format(start, "X"*60, stop)
                self.run_both_dirs(1, seq)

    def test_single_contained(self):
        self.run_both_dirs(1, "ATG"+"X"*60+"TAG")
        self.run_both_dirs(1, "XXXATG"+"X"*60+"TAG")
        self.run_both_dirs(1, "ATG"+"X"*60+"TAGXXX")
        self.run_both_dirs(1, "XXXATG"+"X"*60+"TAGXXX")

    def test_start_without_end(self):
        self.run_both_dirs(1, "XXXATGXXX")

    def test_end_without_start(self):
        self.run_both_dirs(1, "XXXTAGXXX")

    def test_multiple(self):
        # start, stop, start, stop
        self.run_both_dirs(2, "ATG"+"X"*60+"TAGGTG"+"X"*60+"TGA")
        # start, stop, start
        self.run_both_dirs(2, "ATG"+"X"*60+"TAGGTG")
        # stop, start
        self.run_both_dirs(2, "TAGGTGXXX")

    def test_multi_start_single_stop(self):
        self.run_both_dirs(1, "ATGXXXATG"+"X"*60+"TAG")

    def test_single_start_multi_stop(self):
        self.run_both_dirs(1, "ATG"+"X"*60+"TAGXXXTAG")

    def test_interleaved(self):
        self.run_both_dirs(2, "ATGXATGXX"+"X"*60+"TAGXTAG")

class TestOrfLocations(unittest.TestCase):
    def test_contained(self):
        seq = "ATG"+"X"*60+"TAG"
        for offset in [0, 6, 305]:
            if offset == 0:
                orfs = scan_orfs(seq, 1)
            else:
                orfs = scan_orfs(seq, 1, offset=offset)
            assert len(orfs) == 1
            orf = orfs[0]
            assert isinstance(orf.start, ExactPosition)
            assert isinstance(orf.end, ExactPosition)
            assert orf.start == 0 + offset
            assert orf.end == 65 + offset


    def test_start_without_end(self):
        seq = "XXXATGXXX"
        for offset in [0, 6, 305]:
            if offset == 0:
                orfs = scan_orfs(seq, 1)
            else:
                orfs = scan_orfs(seq, 1, offset=offset)
            assert len(orfs) == 1
            orf = orfs[0]
            assert isinstance(orf.start, ExactPosition)
            assert isinstance(orf.end, AfterPosition)
            assert orf.start == 3 + offset

    def test_stop_without_start(self):
        seq = "XXXTAGXXX"
        for offset in [0, 6, 305]:
            if offset == 0:
                orfs = scan_orfs(seq, 1)
            else:
                orfs = scan_orfs(seq, 1, offset=offset)
            assert len(orfs) == 1
            orf = orfs[0]
            assert isinstance(orf.start, BeforePosition)
            assert isinstance(orf.end, ExactPosition)
            assert orf.end == 5 + offset
