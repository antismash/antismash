# License: GNU Affero General Public License v3 or later
# A copy of GNU AGPL v3 should have been included in this software package in LICENSE.txt.

# for test files, silence irrelevant and noisy pylint warnings
# pylint: disable=no-self-use,protected-access,missing-docstring

import unittest

from Bio.Alphabet import generic_dna
from Bio.Seq import Seq
from minimock import mock, restore, TraceTracker, assert_same_trace

import antismash
from antismash.common import path, subprocessing  # mocked, pylint: disable=unused-import
from antismash.common.test.helpers import FakeHSP, DummyRecord, DummyCDS
from antismash.config import build_config
from antismash.detection.full_hmmer import full_hmmer as consts
from antismash.detection import full_hmmer


def _create_dummy_record(reverse=False):
    seq = Seq('GTGGAGCGGTACTAAATGTACTCCACTATCTGCTGATTGGAAACCACGGAGCGCTCTTAG',
              generic_dna)
    strand = 1
    if reverse:
        seq = seq.reverse_complement()
        strand = -1
    rec = DummyRecord(seq=str(seq))

    idx = 1
    for start, end in [(0, 15), (15, 36), (36, 60)]:
        if reverse:
            start, end = len(seq) - end + 3, len(seq) - start  # TODO: check this
        rec.add_cds_feature(DummyCDS(start, end, strand=strand,
                            locus_tag="orf%04d" % idx))
        idx += 1

    return rec


class TestFullhmmer(unittest.TestCase):
    def setUp(self):
        self._old_max_evalue = consts.MAX_EVALUE
        self._old_min_score = consts.MIN_SCORE
        consts.MAX_EVALUE = 0.02
        consts.MIN_SCORE = 1.
        self.config = build_config([], isolated=True,
                                   modules=antismash.get_all_modules())
        self.tracer = TraceTracker()
        mock('antismash.common.path.locate_executable', returns='hmmsearch',
             tracker=self.tracer)
        self.file_list = ['Pfam-A.hmm', 'Pfam-A.hmm.h3f', 'Pfam-A.hmm.h3i',
                          'Pfam-A.hmm.h3m', 'Pfam-A.hmm.h3p']
        mock('antismash.common.path.locate_file', returns_iter=self.file_list,
             tracker=self.tracer)
        mock('antismash.common.subprocessing.run_hmmscan', returns=[])

    def tearDown(self):
        consts.MAX_EVALUE = self._old_max_evalue
        consts.MIN_SCORE = self._old_min_score
        restore()

    def test_check_prereqs(self):
        "Test fullhmmer.check_prereqs()"
        trace = """Called antismash.common.path.locate_executable('hmmscan')
Called antismash.common.path.locate_file(
    '{0}/pfam/Pfam-A.hmm')
Called antismash.common.path.locate_file(
    '{0}/pfam/Pfam-A.hmm.h3f')
Called antismash.common.path.locate_file(
    '{0}/pfam/Pfam-A.hmm.h3i')
Called antismash.common.path.locate_file(
    '{0}/pfam/Pfam-A.hmm.h3m')
Called antismash.common.path.locate_file(
    '{0}/pfam/Pfam-A.hmm.h3p')""".format(self.config.database_dir)
        expected = []
        returned = full_hmmer.check_prereqs()

        self.assertListEqual(expected, returned)
        assert_same_trace(self.tracer, trace)

    def test_check_prereqs_missing_exe(self):
        "Test fullhmmer.check_prereqs() with a missing executable"
        path.locate_executable.mock_returns = None
        trace = """Called antismash.common.path.locate_executable('hmmscan')
Called antismash.common.path.locate_file(
    '{0}/pfam/Pfam-A.hmm')
Called antismash.common.path.locate_file(
    '{0}/pfam/Pfam-A.hmm.h3f')
Called antismash.common.path.locate_file(
    '{0}/pfam/Pfam-A.hmm.h3i')
Called antismash.common.path.locate_file(
    '{0}/pfam/Pfam-A.hmm.h3m')
Called antismash.common.path.locate_file(
    '{0}/pfam/Pfam-A.hmm.h3p')""".format(self.config.database_dir)
        expected = ["Failed to locate file: 'hmmscan'"]
        returned = full_hmmer.check_prereqs()

        self.assertListEqual(expected, returned)
        assert_same_trace(self.tracer, trace)

    def test_check_prereqs_missing_file(self):
        "Test fullhmmer.check_prereqs() with a missing file"
        self.file_list[0] = None
        path.locate_file.mock_returns_iter = self.file_list
        trace = """Called antismash.common.path.locate_executable('hmmscan')
Called antismash.common.path.locate_file(
    '{0}/pfam/Pfam-A.hmm')
Called antismash.common.path.locate_file(
    '{0}/pfam/Pfam-A.hmm.h3f')
Called antismash.common.path.locate_file(
    '{0}/pfam/Pfam-A.hmm.h3i')
Called antismash.common.path.locate_file(
    '{0}/pfam/Pfam-A.hmm.h3m')
Called antismash.common.path.locate_file(
    '{0}/pfam/Pfam-A.hmm.h3p')""".format(self.config.database_dir)
        expected = ["Failed to locate file: 'Pfam-A.hmm'"]
        returned = full_hmmer.check_prereqs()

        self.assertListEqual(expected, returned)
        assert_same_trace(self.tracer, trace)

    def test_calculate_start_and_end(self):
        "Test fullhmmer._calculate_start_end()"
        hsp = FakeHSP(1, 5, 1.5)
        rec = _create_dummy_record()
        feature = rec.get_cds_features()[1]
        start, end = full_hmmer.full_hmmer.calculate_start_and_end(feature, hsp)
        self.assertEqual((18, 30), (start, end))
        self.assertEqual("YSTI", str(rec.seq[start:end].translate()))

        rev = _create_dummy_record(reverse=True)
        feature = rev.get_cds_features()[1]
        start, end = full_hmmer.full_hmmer.calculate_start_and_end(feature, hsp)
        assert start == 30
        assert end == 42
        assert str(rev.seq[start:end].reverse_complement().translate()) == "YSTI"
