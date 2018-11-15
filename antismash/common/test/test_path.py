# License: GNU Affero General Public License v3 or later
# A copy of GNU AGPL v3 should have been included in this software package in LICENSE.txt.

# for test files, silence irrelevant and noisy pylint warnings
# pylint: disable=no-self-use,protected-access,missing-docstring,too-many-public-methods

import os
from tempfile import TemporaryDirectory
import unittest

from antismash.common import path


class TestChangedDir(unittest.TestCase):
    def setUp(self):
        self.original_dir = os.getcwd()

    def tearDown(self):
        os.chdir(self.original_dir)

    def test_functions(self):
        with TemporaryDirectory() as tmpdir:
            assert os.getcwd() != tmpdir
            with path.changed_directory(tmpdir):
                assert os.getcwd() == tmpdir
            assert os.getcwd() == self.original_dir

    def test_with_errors(self):
        try:
            with TemporaryDirectory() as tmpdir:
                with path.changed_directory(tmpdir):
                    raise ValueError("die mid-context")
        except ValueError as err:
            assert "die mid-context" in str(err)
        finally:
            assert os.getcwd() == self.original_dir
