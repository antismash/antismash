# License: GNU Affero General Public License v3 or later
# A copy of GNU AGPL v3 should have been included in this software package in LICENSE.txt.

# for test files, silence irrelevant and noisy pylint warnings
# pylint: disable=no-self-use,protected-access,missing-docstring

import unittest
from unittest.mock import patch

import antismash
from antismash.common import path, subprocessing
from antismash.config import build_config
from antismash.detection import tigrfam


EXPECTED_FILES = ['TIGRFam.hmm', 'TIGRFam.hmm.h3f', 'TIGRFam.hmm.h3i',
                  'TIGRFam.hmm.h3m', 'TIGRFam.hmm.h3p']


@patch.object(subprocessing, "run_hmmscan", return_value=[])
@patch.object(path, "locate_file", side_effect=EXPECTED_FILES)
class TestTIGRFam(unittest.TestCase):
    def setUp(self):
        self._old_max_evalue = tigrfam.MAX_EVALUE
        self._old_min_score = tigrfam.MIN_SCORE
        tigrfam.MAX_EVALUE = 0.02
        tigrfam.MIN_SCORE = 1.
        self.config = build_config([], isolated=True,
                                   modules=antismash.get_all_modules())

    def tearDown(self):
        tigrfam.MAX_EVALUE = self._old_max_evalue
        tigrfam.MIN_SCORE = self._old_min_score

    def test_check_prereqs(self, _patched_locate, _patched_run):
        expected = []
        returned = tigrfam.check_prereqs(self.config)

        self.assertListEqual(expected, returned)

    def test_check_prereqs_missing_exe(self, _patched_locate, _patched_run):
        self.config.executables.__dict__.pop("hmmscan")
        expected = ["Failed to locate executable: 'hmmscan'"]
        returned = tigrfam.check_prereqs(self.config)

        self.assertListEqual(expected, returned)

    def test_check_prereqs_missing_file(self, patched_locate, _patched_run):
        patched_locate.side_effect = [None] + EXPECTED_FILES[1:]
        expected = [f"Failed to locate TIGRFam db in {self.config.database_dir}/tigrfam"]
        returned = tigrfam.check_prereqs(self.config)

        self.assertListEqual(expected, returned)
