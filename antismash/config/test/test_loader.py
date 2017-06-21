# License: GNU Affero General Public License v3 or later
# A copy of GNU AGPL v3 should have been included in this software package in LICENSE.txt.

import os
import unittest

from antismash.config import loader
from antismash.config.args import Config

class TestLoader(unittest.TestCase):
    def setUp(self):
        loader._basedir = os.path.join(os.path.dirname(__file__), 'data')
        loader._user_file_name = os.devnull
        loader._instance_file_name = 'dummy_instance.cfg'

    def tearDown(self):
        Config().__dict__.clear()

    def test_simple(self):
        assert len(Config().__dict__) == 0
        loader.update_config_from_file()
        assert len(Config().__dict__) == 2
        assert Config().base == 'base_value'
        assert Config().child.base == 'child_base_value'
        assert Config().child.bool == False

    def test_override(self):
        Config({'base' : 'other'})
        assert Config().base == 'other'
        loader.update_config_from_file()
        assert Config().base == 'base_value'

    def test_retaining(self):
        Config({'other' : 'value'})
        assert Config().other == 'value'
        loader.update_config_from_file()
        assert Config().other == 'value'
