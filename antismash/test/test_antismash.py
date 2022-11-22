# License: GNU Affero General Public License v3 or later
# A copy of GNU AGPL v3 should have been included in this software package in LICENSE.txt.

# for test files, silence irrelevant and noisy pylint warnings
# pylint: disable=use-implicit-booleaness-not-comparison,protected-access,missing-docstring

from argparse import Namespace
import logging
import os
import unittest
from unittest.mock import call, patch

from antismash import main
from antismash.common.test.helpers import DummyRecord
from antismash.config import build_config, destroy_config, get_config
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
        for option in ["--list-plugins"]:
            options = build_config([option], isolated=False, modules=self.all_modules)
            # don't bother with executables, that's not relevant for this
            with patch.object(main, "_log_found_executables", return_value=None):
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

    def test_canonical_base_filename(self):
        options = build_parser(modules=self.all_modules).parse_args([])
        expected = os.path.join("out", "foo.1_example")
        res = main.canonical_base_filename("foo.1_example.gbk", "out", options)
        assert res == expected
        assert get_config().output_basename == os.path.basename(expected)

        res = main.canonical_base_filename("/some/long/path/foo.1_example.gbff", "out", options)
        assert res == expected

        res = main.canonical_base_filename("foo.1_example.fa", "out", options)
        assert res == expected

        res = main.canonical_base_filename("foo.1_example.gbff.gz", "out", options)
        assert res == expected

        options = build_parser(modules=self.all_modules).parse_args(["--output-basename", "foo.1"])
        expected = os.path.join("out", "foo.1")
        res = main.canonical_base_filename("foo.1_example.gbk", "out", options)
        assert res == expected

        res = main.canonical_base_filename("foo.1_example.gbff", "out", options)
        assert res == expected

        res = main.canonical_base_filename("foo.1_example.fa", "out", options)
        assert res == expected

        res = main.canonical_base_filename("foo.1_example.gbff.gz", "out", options)
        assert res == expected


@patch.object(logging, 'debug')
class TestTimingsLog(unittest.TestCase):
    def test_multi_record(self, mocked_logging):
        results = {"r1": {"a": 2, "b": 4}, "r2": {"a": 2}, "r3": {"b": 3}}
        main.log_module_runtimes(results)
        assert mocked_logging.call_count == 3
        expected = [
            call('Total times taken by modules',),
            call('  %s: %.1fs', 'a', 4.0),
            call('  %s: %.1fs', 'b', 7.0),
        ]
        assert mocked_logging.call_args_list == expected

    def test_nothing_to_report(self, mocked_logging):
        results = {"r1": {}}
        main.log_module_runtimes(results)
        mocked_logging.assert_not_called()
