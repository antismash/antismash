# License: GNU Affero General Public License v3 or later
# A copy of GNU AGPL v3 should have been included in this software package in LICENSE.txt.

# for test files, silence irrelevant and noisy pylint warnings
# pylint: disable=no-self-use,protected-access,missing-docstring

import os
import unittest

from antismash.config import loader
from antismash.config import destroy_config


class TestLoader(unittest.TestCase):
    def setUp(self):
        self.orig_basedir = loader._BASEDIR
        loader._BASEDIR = os.path.join(os.path.dirname(__file__), 'data')

    def tearDown(self):
        loader._BASEDIR = self.orig_basedir
        destroy_config()

    def test_simple(self):
        print(loader.load_config_from_file().__dict__)
        assert loader.load_config_from_file().base == "base_value"
        assert loader.load_config_from_file().child.base == "child_base_value"
