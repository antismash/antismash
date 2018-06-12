# License: GNU Affero General Public License v3 or later
# A copy of GNU AGPL v3 should have been included in this software package in LICENSE.txt.

# for test files, silence irrelevant and noisy pylint warnings
# pylint: disable=no-self-use,protected-access,missing-docstring

import unittest

from Bio.SeqFeature import BeforePosition, AfterPosition, ExactPosition, FeatureLocation

from antismash.common.all_orfs import find_all_orfs, scan_orfs

from .helpers import DummyRecord


class TestOrfCounts(unittest.TestCase):
    def run_both_dirs(self, expected, seq):
        def reverse_location(location, length):
            return location._flip(length)

        record = DummyRecord(seq=seq)
        assert expected == [feat.location for feat in find_all_orfs(record)]
        record.seq = record.seq.reverse_complement()

        expected = [reverse_location(loc, len(seq)) for loc in expected[::-1]]
        assert expected == [feat.location for feat in find_all_orfs(record)]

    def test_empty_sequence(self):
        self.run_both_dirs([], "")

    def test_no_hits(self):
        # no starts or stops
        self.run_both_dirs([], "NNN")

        # nothing > 60
        self.run_both_dirs([], "ATGNNNTGA")
        self.run_both_dirs([], "ATG" + "N"*54 + "TGA")

    def test_all_combos(self):
        expected = [FeatureLocation(ExactPosition(0), ExactPosition(66), strand=1)]
        for start in ('ATG', 'GTG', 'TTG'):
            for stop in ('TAA', 'TAG', 'TGA'):
                seq = "{}{}{}".format(start, "N"*60, stop)
                self.run_both_dirs(expected, seq)

    def test_single_contained(self):
        expected = [FeatureLocation(ExactPosition(0), ExactPosition(66), strand=1)]
        self.run_both_dirs(expected, "ATG"+"N"*60+"TAG")
        self.run_both_dirs(expected, "ATG"+"N"*60+"TAGNNN")
        expected = [FeatureLocation(ExactPosition(3), ExactPosition(69), strand=1)]
        self.run_both_dirs(expected, "NNNATG"+"N"*60+"TAG")
        self.run_both_dirs(expected, "NNNATG"+"N"*60+"TAGNNN")

    def test_start_without_end(self):
        expected = [FeatureLocation(ExactPosition(3), AfterPosition(9), strand=1)]
        self.run_both_dirs(expected, "NNNATGNNN")

    def test_end_without_start(self):
        expected = [FeatureLocation(BeforePosition(0), ExactPosition(6), strand=1)]
        self.run_both_dirs(expected, "NNNTAGNNN")

    def test_multiple(self):
        # start, stop, start, stop
        expected = [FeatureLocation(ExactPosition(0), ExactPosition(66), strand=1),
                    FeatureLocation(ExactPosition(66), ExactPosition(132), strand=1)]
        self.run_both_dirs(expected, "ATG"+"N"*60+"TAGGTG"+"N"*60+"TGA")
        # start, stop, start
        expected[1] = FeatureLocation(ExactPosition(66), AfterPosition(69), strand=1)
        self.run_both_dirs(expected, "ATG"+"N"*60+"TAGGTG")
        # stop, start
        expected = [FeatureLocation(BeforePosition(0), ExactPosition(3), strand=1),
                    FeatureLocation(ExactPosition(3), AfterPosition(9), strand=1)]
        self.run_both_dirs(expected, "TAGGTGNNN")

    def test_multi_start_single_stop(self):
        seq = "ATGNNNATG" + "N"*60 + "TAG"
        expected = [FeatureLocation(ExactPosition(0), ExactPosition(72), strand=1)]
        assert expected == [feat.location for feat in find_all_orfs(DummyRecord(seq=seq))]
        seq = str(DummyRecord(seq=seq).seq.reverse_complement())
        expected[0].strand = -1
        assert expected == [feat.location for feat in find_all_orfs(DummyRecord(seq=seq))]

    def test_single_start_multi_stop(self):
        seq = "ATG"+"N"*60+"TAGNNNTAG"
        expected = [FeatureLocation(ExactPosition(0), ExactPosition(66), strand=1)]
        assert expected == [feat.location for feat in find_all_orfs(DummyRecord(seq=seq))]
        seq = str(DummyRecord(seq=seq).seq.reverse_complement())
        expected = [FeatureLocation(ExactPosition(6), ExactPosition(72), strand=-1)]
        assert expected == [feat.location for feat in find_all_orfs(DummyRecord(seq=seq))]

    def test_interleaved(self):
        expected = [FeatureLocation(ExactPosition(0), ExactPosition(72), strand=1),
                    FeatureLocation(ExactPosition(4), ExactPosition(76), strand=1)]
        self.run_both_dirs(expected, "ATGNATGNN"+"N"*60+"TAGNTAG")

    def test_multiframe(self):
        starts = [3, 25, 72, 77]
        seq = ["X"] * 200
        for start in starts:
            seq[start:start + 3] = "ATG"
            seq[start + 63: start + 66] = "TAA"
        assert len(seq) == 200
        result = scan_orfs("".join(seq), direction=1)
        assert sorted([orf.start for orf in result]) == starts
        assert len(result) == 4

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
            assert orf.end == 66 + offset

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
            assert orf.end == 6 + offset
