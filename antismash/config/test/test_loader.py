# License: GNU Affero General Public License v3 or later
# A copy of GNU AGPL v3 should have been included in this software package in LICENSE.txt.

import os
import unittest

from antismash.config import loader
from antismash.config.args import Config

class TestLoader(unittest.TestCase):
    def setUp(self):
# pylint: disable=protected-access
        self.orig_basedir = loader._BASEDIR
        self.orig_user = loader._USER_FILE_NAME
        self.orig_instance = loader._INSTANCE_FILE_NAME
        loader._BASEDIR = os.path.join(os.path.dirname(__file__), 'data')
        loader._USER_FILE_NAME = os.devnull
        loader._INSTANCE_FILE_NAME = 'dummy_instance.cfg'
# pylint: enable=protected-access

    def tearDown(self):
# pylint: disable=protected-access
        loader._BASEDIR = self.orig_basedir
        loader._USER_FILE_NAME = self.orig_user
        loader._INSTANCE_FILE_NAME = self.orig_instance
        Config().__dict__.clear()
# pylint: enable=protected-access

    def test_simple(self):
        assert not Config()
        loader.update_config_from_file()
        assert len(Config()) == 2
        assert Config().base == 'base_value'
        assert Config().child.base == 'child_base_value'
        assert Config().child.bool is False

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
