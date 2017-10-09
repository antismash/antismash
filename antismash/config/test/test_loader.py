# License: GNU Affero General Public License v3 or later
# A copy of GNU AGPL v3 should have been included in this software package in LICENSE.txt.

# for test files, silence irrelevant and noisy pylint warnings
# pylint: disable=no-self-use,protected-access,missing-docstring

import os
import unittest

from antismash.config import loader
from antismash.config import update_config, get_config, destroy_config


class TestLoader(unittest.TestCase):
    def setUp(self):
        self.orig_basedir = loader._BASEDIR
        self.orig_user = loader._USER_FILE_NAME
        self.orig_instance = loader._INSTANCE_FILE_NAME
        loader._BASEDIR = os.path.join(os.path.dirname(__file__), 'data')
        loader._USER_FILE_NAME = os.devnull
        loader._INSTANCE_FILE_NAME = 'dummy_instance.cfg'

    def tearDown(self):
        loader._BASEDIR = self.orig_basedir
        loader._USER_FILE_NAME = self.orig_user
        loader._INSTANCE_FILE_NAME = self.orig_instance
        destroy_config()

    def test_simple(self):
        assert not get_config()
        loader.update_config_from_file()
        assert len(get_config()) == 2
        assert get_config().base == 'base_value'
        assert get_config().child.base == 'child_base_value'
        assert get_config().child.bool is False

    def test_override(self):
        update_config({'base': 'other'})
        assert get_config().base == 'other'
        loader.update_config_from_file()
        assert get_config().base == 'base_value'

    def test_retaining(self):
        update_config({'other': 'value'})
        assert get_config().other == 'value'
        loader.update_config_from_file()
        assert get_config().other == 'value'
