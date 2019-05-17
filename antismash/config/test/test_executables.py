# License: GNU Affero General Public License v3 or later
# A copy of GNU AGPL v3 should have been included in this software package in LICENSE.txt.

# for test files, silence irrelevant and noisy pylint warnings
# pylint: disable=no-self-use,protected-access,missing-docstring

from collections import OrderedDict
import os
import unittest
from unittest.mock import patch

from antismash.config import executables


class TestGetPaths(unittest.TestCase):
    def test_bad_formats(self):
        for bad in ["a", "a=", "=a", "a=a=", "a=a=a"]:
            with self.assertRaisesRegex(ValueError, "invalid .* format"):
                executables.get_executable_paths(bad)

    @patch.object(executables, "find_executable_path", return_value="/some/other/abs/path")
    @patch.object(os.path, "exists", return_value=True)
    def test_duplicated(self, _mocked_find, _mocked_exists):
        paths = ["diamond=some_alias", "diamond=/some/abs/path"]
        for direction in [1, -1]:
            arg = ",".join(paths[::direction])
            with self.assertRaisesRegex(ValueError, "multiple paths specified"):
                executables.get_executable_paths(arg)

    @patch.object(executables, "find_executable_path", return_value="/abs/path/alt_diamond")
    @patch("os.path.exists", return_value=True)
    def test_duplicated_but_same_value(self, _mocked_find, _mocked_exists):
        paths = ["diamond=alt_diamond", "diamond=/abs/path/alt_diamond"]
        for direction in [1, -1]:
            arg = ",".join(paths[::direction])
            assert executables.get_executable_paths(arg)["diamond"] == "/abs/path/alt_diamond"

    @patch("os.path.exists", return_value=False)
    def test_missing_abs_target(self, _mocked_exists):
        with self.assertRaisesRegex(ValueError, "no such file"):
            executables.get_executable_paths("diamond=/bad/path")

    @patch.object(executables, "find_executable_path", return_value=None)
    def test_unlocatable(self, _mocked_find):
        with self.assertRaisesRegex(ValueError, "cannot find"):
            executables.get_executable_paths("diamond=wrong_alias")

    @patch.object(executables, "_ALTERNATE_EXECUTABLE_NAMES",
                  OrderedDict([("missing", ["missing"]), ("found", ["found"])]))
    @patch.object(executables, "find_executable_path", side_effect=["", "/found"])
    def test_missing_default_path(self, _mocked_find):
        defaults = executables.get_default_paths()
        assert "missing" not in defaults
        assert defaults["found"] == "/found"
