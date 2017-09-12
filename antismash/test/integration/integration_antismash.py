# License: GNU Affero General Public License v3 or later
# A copy of GNU AGPL v3 should have been included in this software package in LICENSE.txt.

import os
import sys
import unittest

from antismash.main import run_antismash, get_all_modules
from antismash.config.args import build_parser, Config

class TestAntismash(unittest.TestCase):
    def setUp(self):
        args = ["run_antismash.py"]
        self.parser = build_parser(modules=get_all_modules())
        self.default_options = self.parser.parse_args(args)
        self.default_options.minimal = True
        self.default_options.all_enabled_modules = []
        Config(self.default_options)

    def tearDown(self):
        Config().__dict__.clear()

    def test_nisin_minimal(self):
        path = os.path.join(os.path.dirname(__file__), "data", "nisin.gbk")
        assert Config().minimal
        run_antismash(path, Config())
