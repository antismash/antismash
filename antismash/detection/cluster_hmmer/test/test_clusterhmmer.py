# License: GNU Affero General Public License v3 or later
# A copy of GNU AGPL v3 should have been included in this software package in LICENSE.txt.

# for test files, silence irrelevant and noisy pylint warnings
# pylint: disable=use-implicit-booleaness-not-comparison,protected-access,missing-docstring

import unittest
from unittest.mock import patch

import antismash
from antismash.common import path, pfamdb, subprocessing
from antismash.config import build_config
from antismash.detection import cluster_hmmer


EXPECTED_FILES = ['Pfam-A.hmm', 'Pfam-A.hmm.h3f', 'Pfam-A.hmm.h3i',
                  'Pfam-A.hmm.h3m', 'Pfam-A.hmm.h3p']


@patch.object(subprocessing, "run_hmmscan", return_value=[])
@patch.object(path, "locate_file", side_effect=EXPECTED_FILES)
class TestClusterhmmer(unittest.TestCase):
    def setUp(self):
        self._old_max_evalue = cluster_hmmer.MAX_EVALUE
        self._old_min_score = cluster_hmmer.MIN_SCORE
        cluster_hmmer.MAX_EVALUE = 0.02
        cluster_hmmer.MIN_SCORE = 1.
        self.config = build_config([], isolated=True,
                                   modules=antismash.get_all_modules())
        self.latest_pfam = pfamdb.find_latest_database_version(self.config.database_dir)

    def tearDown(self):
        cluster_hmmer.MAX_EVALUE = self._old_max_evalue
        cluster_hmmer.MIN_SCORE = self._old_min_score

    def test_check_prereqs(self, _patched_locate, _patched_run):
        expected = []
        returned = cluster_hmmer.check_prereqs(self.config)

        self.assertListEqual(expected, returned)

    def test_check_prereqs_missing_exe(self, _patched_locate, _patched_run):
        self.config.executables.__dict__.pop("hmmscan")
        expected = ["Failed to locate executable: 'hmmscan'"]
        returned = cluster_hmmer.check_prereqs(self.config)

        self.assertListEqual(expected, returned)

    def test_check_prereqs_missing_file(self, patched_locate, _patched_run):
        patched_locate.side_effect = [None] + EXPECTED_FILES[1:]
        expected = [f"Failed to locate file: 'Pfam-A.hmm' in {self.config.database_dir}/pfam/{self.latest_pfam}"]
        returned = cluster_hmmer.check_prereqs(self.config)

        self.assertListEqual(expected, returned)
