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
    mock_run_diamond.assert_called_once_with("version", use_default_opts=False)


@patch.object(diamond, "run_diamond", return_value=DummyResult("lots of useless text"))
@patch.object(diamond, "run_diamond_version", return_value="diamond version X.Y.Z")
def test_diamond_makedb(_mocked_version, mock_run_diamond):
    subprocessing.run_diamond_makedb("fake.dmnd", "fake.fasta")
    assert mock_run_diamond.call_count == 1
    assert mock_run_diamond.call_args.kwargs == {"use_default_opts": True}
    assert mock_run_diamond.call_args[0] == ("makedb", ["--db", "fake.dmnd", "--in", "fake.fasta"])


@patch.object(diamond, "TemporaryDirectory")
def test_diamond_tmpdir_bug(mocked_tempdir):
    temp_dir = "some_temp_dir"
    mocked_tempdir.return_value.__enter__.return_value = temp_dir
    # make sure the mock works before continuing, since it's a bit complicated
    with diamond.TemporaryDirectory() as temp:
        assert temp == temp_dir

    def run(version):
        with patch.object(diamond, "execute") as mocked_execute:
            with patch.object(diamond, "run_diamond_version", return_value=version):
                subprocessing.run_diamond_makedb("fake.dmnd", "fake.fasta")
                assert mocked_execute.call_count == 1
                args = mocked_execute.call_args[0][0]
        return args

    # known good boundaries and some odd formats that can't be predicted
    for version in ["0.9.3", "2.0.4", "2.1.7", "vX.Y.Z"]:
        args_used = run(version)
        assert "--tmpdir" in args_used
    # known bad boundaries
    for version in ["2.1.0", "2.1.4", "2.1.6"]:
        args_used = run(version)
        assert "--tmpdir" not in args_used


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


def test_output_parsing():
    expected = [
        diamond.Hit(query_id="queryid", reference_id="refid1", identity=95, evalue=1e-20, bitscore=250.5),
        diamond.Hit(query_id="queryid", reference_id="refid2", identity=85, evalue=1e-25, bitscore=151.2),
    ]

    content = (
        "queryid\trefid1\t95\t20\t1\t2\t80\t100\t20\t40\t1e-20\t250.5\n"
        "queryid\trefid2\t85\t20\t3\t2\t5\t25\t25\t45\t1e-25\t151.2\n"
    )
    hits = diamond.parse_tabular_output(content)
    assert hits == expected
