# License: GNU Affero General Public License v3 or later
# A copy of GNU AGPL v3 should have been included in this software package in LICENSE.txt.

# for test files, silence irrelevant and noisy pylint warnings
# pylint: disable=no-self-use,protected-access,missing-docstring

import unittest

from Bio.SeqFeature import ExactPosition, FeatureLocation

from antismash.common.all_orfs import (
    find_all_orfs,
    find_intergenic_areas,
    get_trimmed_orf,
    scan_orfs,
)

from .helpers import DummyCDS, DummyRecord


class TestOrfCounts(unittest.TestCase):
    def run_both_dirs(self, expected, seq, min_length=60):
        def reverse_location(location, length):
            return location._flip(length)

        record = DummyRecord(seq=seq)
        assert expected == [feat.location for feat in find_all_orfs(record, min_length=min_length)]
        record.seq = record.seq.reverse_complement()
        print(record.seq)

        expected = [reverse_location(loc, len(seq)) for loc in expected[::-1]]
        assert expected == [feat.location for feat in find_all_orfs(record, min_length=min_length)]

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
        self.run_both_dirs([], "NNNATGNNN")

    def test_end_without_start(self):
        self.run_both_dirs([], "NNNTAGNNN")

    def test_multiple(self):
        # start, stop, start, stop
        expected = [FeatureLocation(ExactPosition(0), ExactPosition(66), strand=1),
                    FeatureLocation(ExactPosition(66), ExactPosition(132), strand=1)]
        self.run_both_dirs(expected, "ATG"+"N"*60+"TAGGTG"+"N"*60+"TGA")
        # start, stop, start
        expected.pop()
        self.run_both_dirs(expected, "ATG"+"N"*60+"TAGGTG")
        # stop, start
        expected = []
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


class TestIntegenic(unittest.TestCase):
    def test_no_cdses(self):
        assert find_intergenic_areas(0, 120, []) == [(0, 120)]
        assert find_intergenic_areas(20, 100, []) == [(20, 100)]

    def test_min_length(self):
        cdses = [DummyCDS(30, 36, strand=1), DummyCDS(39, 45, strand=-1)]
        assert find_intergenic_areas(0, 60, cdses, min_length=6) == [(0, 30), (45, 60)]
        assert find_intergenic_areas(0, 60, cdses, min_length=3) == [(0, 30), (36, 39), (45, 60)]
        assert find_intergenic_areas(0, 60, cdses, min_length=20) == [(0, 30)]

    def test_simple(self):
        cdses = [DummyCDS(30, 36, strand=1), DummyCDS(39, 45, strand=-1)]
        areas = find_intergenic_areas(0, 120, []) == [(0, 30), (36, 39), (45, 120)]

    def test_padding(self):
        cdses = [DummyCDS(30, 45, strand=1)]
        areas = find_intergenic_areas(0, 60, cdses, padding=0) == [(0, 30), (45, 60)]
        areas = find_intergenic_areas(0, 60, cdses, padding=3) == [(0, 27), (48, 60)]
        areas = find_intergenic_areas(12, 50, cdses, padding=3) == [(12, 27), (48, 50)]

    def test_overlapping_cds(self):
        cdses = [DummyCDS(30, 42, strand=1), DummyCDS(41, 50, strand=-1)]
        areas = find_intergenic_areas(0, 60, cdses) == [(0, 30), (50, 60)]


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
            assert not orfs

    def test_stop_without_start(self):
        seq = "NNNTAGNNN"
        for offset in [0, 6, 305]:
            if offset == 0:
                orfs = scan_orfs(seq, 1)
            else:
                orfs = scan_orfs(seq, 1, offset=offset)
            assert not orfs


class TestOrfTrimming(unittest.TestCase):
    def setUp(self):
        self.record = DummyRecord(seq="nnnATGnnnATGnnnATGnnnTGAnnnnnn")
        self.cds = DummyCDS(3, 24, strand=1)

    def test_bad_args(self):
        with self.assertRaisesRegex(ValueError, "minimum .* cannot be greater than maximum"):
            get_trimmed_orf(self.cds, self.record, min_length=100, max_length=99)

    def test_simple(self):
        new = get_trimmed_orf(self.cds, self.record)
        assert new
        assert new.location.start == 15
        assert new.location.start >= self.cds.location.start
        assert new.location.end == self.cds.location.end
        assert new.location.strand == self.cds.location.strand
        assert new.translation == "MX"

    def test_includes(self):
        for include, expected in [(4, 3), (13, 15)]:  # coordinate _within_ the CDS, new within record
            new = get_trimmed_orf(self.cds, self.record, include=include)
            assert new
            assert new.location.start >= self.cds.location.start
            start = new.location.start
            end = new.location.end
            print(include, f"{self.record.seq[:start]}{self.record.seq[start:end].upper()}{self.record.seq[end:]}")
            print(new.location.extract(self.record.seq))
            assert new.location.start == expected
            assert new.location.end == self.cds.location.end
            assert new.location.strand == self.cds.location.strand

    def test_reverse_strand(self):
        record = DummyRecord(seq=self.record.seq.reverse_complement())
        cds = DummyCDS(6, 27, strand=-1)
        new = get_trimmed_orf(cds, record)
        assert new.location.start == cds.location.start
        assert new.location.end == 15
        assert new.location.strand == cds.location.strand

    def test_reverse_include(self):
        record = DummyRecord(seq=self.record.seq.reverse_complement())
        cds = DummyCDS(6, 27, strand=-1)
        new = get_trimmed_orf(cds, record, include=10)
        assert new.location.start == cds.location.start
        assert new.location.end == 21
        assert new.location.strand == cds.location.strand

    def test_label(self):
        new = get_trimmed_orf(self.cds, self.record)
        assert new
        assert new.locus_tag.startswith("allorf")
        new = get_trimmed_orf(self.cds, self.record, label="newlabel")
        assert new.locus_tag == "newlabel"

    def test_not_possible(self):
        # no start codon within 6 bases of stop
        assert get_trimmed_orf(self.cds, self.record, max_length=6) is None

        # min length longer than CDS (i.e. no start codons available)
        assert get_trimmed_orf(self.cds, self.record, min_length=len(self.cds.location) + 1,
                               max_length=100) is None

        # max length falls before the include coordinate
        assert get_trimmed_orf(self.cds, self.record, include=0, max_length=6) is None
