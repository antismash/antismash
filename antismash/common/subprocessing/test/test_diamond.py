# License: GNU Affero General Public License v3 or later
# A copy of GNU AGPL v3 should have been included in this software package in LICENSE.txt.

# for test files, silence irrelevant and noisy pylint warnings
# pylint: disable=use-implicit-booleaness-not-comparison,protected-access,missing-docstring

import unittest
from unittest.mock import patch

from antismash.common import path, subprocessing
from antismash.common.subprocessing import diamond

from .helpers import DummyResult


@patch("antismash.common.subprocessing.diamond.run_diamond", return_value=DummyResult("diamond version 1.2.3"))
def test_diamond_version(mock_run_diamond):
    version = subprocessing.run_diamond_version()
    assert version == "1.2.3"
    mock_run_diamond.assert_called_once_with("version")


@patch("antismash.common.subprocessing.diamond.run_diamond", return_value=DummyResult("lots of useless text"))
def test_diamond_makedb(mock_run_diamond):
    subprocessing.run_diamond_makedb("fake.dmnd", "fake.fasta")
    mock_run_diamond.assert_called_once_with("makedb", ["--db", "fake.dmnd", "--in", "fake.fasta"])


class TestDiamondDatabaseChecks(unittest.TestCase):
    def setUp(self):
        self.format0_file = path.get_full_path(__file__, "data", "format0.dmnd")
        self.format1_file = path.get_full_path(__file__, "data", "format1.dmnd")
        self.empty = path.get_full_path(__file__, "data", "empty.dmnd")

    def test_extract_db_format(self):
        assert diamond._extract_db_format(self.format0_file) == 0
        assert diamond._extract_db_format(self.format1_file) == 1
        with self.assertRaises(ValueError):
            diamond._extract_db_format(self.empty)
