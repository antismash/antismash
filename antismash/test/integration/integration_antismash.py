# License: GNU Affero General Public License v3 or later
# A copy of GNU AGPL v3 should have been included in this software package in LICENSE.txt.

# for test files, silence irrelevant and noisy pylint warnings
# pylint: disable=no-self-use,protected-access,missing-docstring

import os
from tempfile import TemporaryDirectory
import unittest

from antismash.main import run_antismash, get_all_modules
from antismash.config import build_config, get_config
from antismash.config.args import build_parser


class TestAntismash(unittest.TestCase):
    def setUp(self):
        args = ["run_antismash.py", "--minimal"]
        self.parser = build_parser(modules=get_all_modules())
        self.default_options = build_config(args, parser=self.parser)
        self.default_options.all_enabled_modules = []
        self.temp_dir = TemporaryDirectory()
        self.default_options.output_dir = self.temp_dir.name

    def tearDown(self):
        get_config().__dict__.clear()
        self.temp_dir.cleanup()

    def test_nisin_minimal(self):
        path = os.path.abspath(os.path.join(os.path.dirname(__file__), "data", "nisin.gbk"))
        run_antismash(path, self.default_options)
