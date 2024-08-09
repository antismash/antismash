# License: GNU Affero General Public License v3 or later
# A copy of GNU AGPL v3 should have been included in this software package in LICENSE.txt.

# for test files, silence irrelevant and noisy pylint warnings
# pylint: disable=use-implicit-booleaness-not-comparison,protected-access,missing-docstring

import os
from tempfile import NamedTemporaryFile
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

    def test_alternate(self):
        assert not loader.load_config_from_file("/dev/null").__dict__
        original_path = os.path.join(loader._BASEDIR, "default.cfg")
        with NamedTemporaryFile() as altered:
            with open(original_path, encoding="utf-8") as original:
                content = original.read().replace("child_base_value", "altered_val")
                altered.write(content.encode())
            altered.seek(0)
            config = loader.load_config_from_file(altered.name)
        assert config.child.base == "altered_val"
