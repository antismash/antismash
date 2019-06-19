# License: GNU Affero General Public License v3 or later
# A copy of GNU AGPL v3 should have been included in this software package in LICENSE.txt.

# for test files, silence irrelevant and noisy pylint warnings
# pylint: disable=no-self-use,protected-access,missing-docstring

from argparse import Namespace
import logging  # pylint: disable=unused-import
import unittest

from minimock import mock, restore, TraceTracker, assert_same_trace

from antismash import main
from antismash.common.test.helpers import DummyRecord
from antismash.config import build_config, destroy_config
from antismash.config.args import build_parser
from antismash.detection import cluster_hmmer, full_hmmer


class TestAntismash(unittest.TestCase):
    def setUp(self):
        self.all_modules = main.get_all_modules()
        self.default_options = build_parser(modules=self.all_modules).parse_args([])

    def tearDown(self):
        destroy_config()

    def test_default_options(self):
        # default options should work with all normal modules
        # barring those using databases not packed with antismash
        options = self.default_options
        modules = list(self.all_modules)
        for special in [full_hmmer, cluster_hmmer]:
            if special in modules:
                modules.pop(modules.index(special))
        assert main.verify_options(options, modules)

    def test_help_options(self):
        for option in ["--list-plugins", "--check-prereqs"]:
            options = build_config([option], isolated=False, modules=self.all_modules)
            ret_val = main.run_antismash("", options)
            assert ret_val == 0

    def test_prepare(self):
        class DummyModule:  # pylint: disable=too-few-public-methods
            def __init__(self):
                self.prep_called = False

            def prepare_data(self):
                self.prep_called = True

        dummy = DummyModule()
        assert not dummy.prep_called
        main.prepare_module_data([dummy])
        assert dummy.prep_called

    def test_antismash_comment(self):
        rec = DummyRecord()
        options = Namespace()
        options.start = -1
        options.end = -1
        options.version = "5.dummy"
        bio = rec.to_biopython()

        main.add_antismash_comments([(rec, bio)], options)
        assert "##antiSMASH-Data-START##" in bio.annotations["comment"]
        assert "##antiSMASH-Data-END##" in bio.annotations["comment"]
        assert "Version" in bio.annotations["comment"] and options.version in bio.annotations["comment"]
        assert "Original ID" not in bio.annotations["comment"]
        assert "Starting at" not in bio.annotations["comment"]
        assert "Ending at" not in bio.annotations["comment"]

        bio.annotations["comment"] = ""
        options.start = 7
        main.add_antismash_comments([(rec, bio)], options)
        assert "Original ID" not in bio.annotations["comment"]
        assert "Starting at  :: 7\n" in bio.annotations["comment"]

        bio.annotations["comment"] = ""
        options.start = -1
        options.end = 1000
        main.add_antismash_comments([(rec, bio)], options)
        assert "Original ID" not in bio.annotations["comment"]
        assert "Ending at    :: 1000\n" in bio.annotations["comment"]

        bio.annotations["comment"] = ""
        options.end = -1
        rec.original_id = "something else"
        main.add_antismash_comments([(rec, bio)], options)
        assert "Original ID" in bio.annotations["comment"] and "something else" in bio.annotations["comment"]


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
