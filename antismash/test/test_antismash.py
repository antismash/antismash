# License: GNU Affero General Public License v3 or later
# A copy of GNU AGPL v3 should have been included in this software package in LICENSE.txt.

# for test files, silence irrelevant and noisy pylint warnings
# pylint: disable=no-self-use,protected-access,missing-docstring

import logging  # pylint: disable=unused-import
import unittest

from minimock import mock, restore, TraceTracker, assert_same_trace

from antismash import main
from antismash.config import build_config, destroy_config
from antismash.config.args import build_parser

class TestAntismash(unittest.TestCase):
    def setUp(self):
        self.all_modules = main.get_all_modules()
        self.default_options = build_parser(modules=self.all_modules).parse_args([])

    def tearDown(self):
        destroy_config()

    def test_default_options(self):
        # default options should work with all normal modules
        options = self.default_options
        assert main.verify_options(options, self.all_modules)

        # adding an incompatibility should not be ok
        options.tta = True
        options.input_type = 'prot'
        assert not main.verify_options(options, self.all_modules)

        # and just to be sure the above isn't just because tta
        options.input_type = 'nucl'
        assert main.verify_options(options, self.all_modules)

    def test_help_options(self):
        for option in ["--list-plugins", "--check-prereqs"]:
            options = build_config([option], isolated=False, modules=self.all_modules)
            ret_val = main.run_antismash("", options)
            assert ret_val == 0

class TestTimingsLog(unittest.TestCase):
    def setUp(self):
        self.trace_tracker = TraceTracker()
        mock('logging.debug', tracker=self.trace_tracker, returns=None)

    def tearDown(self):
        restore()

    def test_multi_record(self):
        results = {"r1": {"a": 2, "b": 4}, "r2": {"a": 2}, "r3": {"b": 3}}
        expected = ("Called logging.debug('Total times taken by modules')\n"
                    "Called logging.debug('  %s: %.1fs', 'a', 4.0)\n"
                    "Called logging.debug('  %s: %.1fs', 'b', 7.0)\n")
        main.log_module_runtimes(results)
        assert_same_trace(self.trace_tracker, expected)

    def test_nothing_to_report(self):
        results = {"r1": {}}
        main.log_module_runtimes(results)
        assert_same_trace(self.trace_tracker, "")
