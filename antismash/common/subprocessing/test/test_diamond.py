# License: GNU Affero General Public License v3 or later
# A copy of GNU AGPL v3 should have been included in this software package in LICENSE.txt.

# for test files, silence irrelevant and noisy pylint warnings
# pylint: disable=no-self-use,protected-access,missing-docstring

from unittest.mock import patch

from antismash.common import subprocessing

class DummyResult(subprocessing.RunResult):
    def __init__(self, stdout: str):
        super().__init__(["dummy"], stdout.encode(), b"", 0, True, True)


@patch("antismash.common.subprocessing.diamond.run_diamond", return_value=DummyResult("diamond version 1.2.3"))
def test_diamond_version(mock_run_diamond):
    version = subprocessing.run_diamond_version()
    assert version == "1.2.3"
    mock_run_diamond.assert_called_once_with("version")


@patch("antismash.common.subprocessing.diamond.run_diamond", return_value=DummyResult("lots of useless text"))
def test_diamond_makedb(mock_run_diamond):
    subprocessing.run_diamond_makedb("fake.dmnd", "fake.fasta")
    mock_run_diamond.assert_called_once_with("makedb",
        ["--db", "fake.dmnd", "--in", "fake.fasta"])
