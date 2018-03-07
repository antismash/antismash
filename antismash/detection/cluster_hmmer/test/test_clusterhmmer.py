# License: GNU Affero General Public License v3 or later
# A copy of GNU AGPL v3 should have been included in this software package in LICENSE.txt.

# for test files, silence irrelevant and noisy pylint warnings
# pylint: disable=no-self-use,protected-access,missing-docstring

import unittest

from minimock import mock, restore, TraceTracker, assert_same_trace

import antismash
from antismash.common import path, pfamdb, subprocessing  # mocked, pylint: disable=unused-import
from antismash.config import build_config
from antismash.detection import cluster_hmmer


class TestFullhmmer(unittest.TestCase):
    def setUp(self):
        self._old_max_evalue = cluster_hmmer.MAX_EVALUE
        self._old_min_score = cluster_hmmer.MIN_SCORE
        cluster_hmmer.MAX_EVALUE = 0.02
        cluster_hmmer.MIN_SCORE = 1.
        self.config = build_config([], isolated=True,
                                   modules=antismash.get_all_modules())
        self.latest_pfam = pfamdb.find_latest_database_version(self.config.database_dir)
        self.tracer = TraceTracker()
        mock('antismash.common.path.locate_executable', returns='hmmsearch',
             tracker=self.tracer)
        self.file_list = ['Pfam-A.hmm', 'Pfam-A.hmm.h3f', 'Pfam-A.hmm.h3i',
                          'Pfam-A.hmm.h3m', 'Pfam-A.hmm.h3p']
        mock('antismash.common.path.locate_file', returns_iter=self.file_list,
             tracker=self.tracer)
        mock('antismash.common.subprocessing.run_hmmscan', returns=[])

        self.expected_trace = """Called antismash.common.path.locate_executable('hmmscan')
Called antismash.common.path.locate_file(
    '{0}/pfam/{1}/Pfam-A.hmm')
Called antismash.common.path.locate_file(
    '{0}/pfam/{1}/Pfam-A.hmm.h3f')
Called antismash.common.path.locate_file(
    '{0}/pfam/{1}/Pfam-A.hmm.h3i')
Called antismash.common.path.locate_file(
    '{0}/pfam/{1}/Pfam-A.hmm.h3m')
Called antismash.common.path.locate_file(
    '{0}/pfam/{1}/Pfam-A.hmm.h3p')""".format(self.config.database_dir, self.latest_pfam)

    def tearDown(self):
        cluster_hmmer.MAX_EVALUE = self._old_max_evalue
        cluster_hmmer.MIN_SCORE = self._old_min_score
        restore()

    def test_check_prereqs(self):
        "Test fullhmmer.check_prereqs()"
        expected = []
        returned = cluster_hmmer.check_prereqs()

        self.assertListEqual(expected, returned)
        assert_same_trace(self.tracer, self.expected_trace)

    def test_check_prereqs_missing_exe(self):
        "Test fullhmmer.check_prereqs() with a missing executable"
        path.locate_executable.mock_returns = None
        expected = ["Failed to locate executable: 'hmmscan'"]
        returned = cluster_hmmer.check_prereqs()

        self.assertListEqual(expected, returned)
        assert_same_trace(self.tracer, self.expected_trace)

    def test_check_prereqs_missing_file(self):
        "Test fullhmmer.check_prereqs() with a missing file"
        self.file_list[0] = None
        path.locate_file.mock_returns_iter = self.file_list
        expected = ["Failed to locate file: 'Pfam-A.hmm' in %s/pfam/%s" % (self.config.database_dir, self.latest_pfam)]
        returned = cluster_hmmer.check_prereqs()

        self.assertListEqual(expected, returned)
        assert_same_trace(self.tracer, self.expected_trace)
