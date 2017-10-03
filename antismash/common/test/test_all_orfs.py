# License: GNU Affero General Public License v3 or later
# A copy of GNU AGPL v3 should have been included in this software package in LICENSE.txt.

# for test files, silence irrelevant and noisy pylint warnings
# pylint: disable=no-self-use,protected-access,missing-docstring

import unittest

from Bio.SeqFeature import BeforePosition, AfterPosition, ExactPosition

from antismash.common.all_orfs import find_all_orfs, scan_orfs

from .helpers import DummyRecord

class TestOrfCounts(unittest.TestCase):
    def run_both_dirs(self, expected, seq):
        assert expected == find_all_orfs(DummyRecord(seq=seq))

        # replace for reverse direction testing
        for start, comp in [('ATG', 'TAC'), ('GTG', 'CAC'), ('TTG', 'AAC')]:
            seq = seq.replace(start, comp)
        for stop, comp in [('TAA', 'ATT'), ('TAG', 'ATC'), ('TGA', 'ACT')]:
            seq = seq.replace(stop, comp)
        assert expected == find_all_orfs(DummyRecord(seq=seq))

    def test_empty_sequence(self):
        self.run_both_dirs(0, "")

    def test_no_hits(self):
        # no starts or stops
        self.run_both_dirs(0, "NNN")

        # nothing > 60
        self.run_both_dirs(0, "ATGNNNTGA")
        self.run_both_dirs(0, "ATG"+ "N"*54 +"TGA")

    def test_all_combos(self):
        for start in ('ATG', 'GTG', 'TTG'):
            for stop in ('TAA', 'TAG', 'TGA'):
                seq = "{}{}{}".format(start, "N"*60, stop)
                self.run_both_dirs(1, seq)

    def test_single_contained(self):
        self.run_both_dirs(1, "ATG"+"N"*60+"TAG")
        self.run_both_dirs(1, "NNNATG"+"N"*60+"TAG")
        self.run_both_dirs(1, "ATG"+"N"*60+"TAGNNN")
        self.run_both_dirs(1, "NNNATG"+"N"*60+"TAGNNN")

    def test_start_without_end(self):
        self.run_both_dirs(1, "NNNATGNNN")

    def test_end_without_start(self):
        self.run_both_dirs(1, "NNNTAGNNN")

    def test_multiple(self):
        # start, stop, start, stop
        self.run_both_dirs(2, "ATG"+"N"*60+"TAGGTG"+"N"*60+"TGA")
        # start, stop, start
        self.run_both_dirs(2, "ATG"+"N"*60+"TAGGTG")
        # stop, start
        self.run_both_dirs(2, "TAGGTGNNN")

    def test_multi_start_single_stop(self):
        self.run_both_dirs(1, "ATGNNNATG"+"N"*60+"TAG")

    def test_single_start_multi_stop(self):
        self.run_both_dirs(1, "ATG"+"N"*60+"TAGNNNTAG")

    def test_interleaved(self):
        self.run_both_dirs(2, "ATGNATGNN"+"N"*60+"TAGNTAG")

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
        seq = "NNNATGNNN"
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
        seq = "NNNTAGNNN"
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
