# License: GNU Affero General Public License v3 or later
# A copy of GNU AGPL v3 should have been included in this software package in LICENSE.txt.

# for test files, silence irrelevant and noisy pylint warnings
# pylint: disable=use-implicit-booleaness-not-comparison,protected-access,missing-docstring

from argparse import Namespace
import logging
import os
from tempfile import NamedTemporaryFile
import unittest
from unittest.mock import call, patch

import pytest

from antismash import __main__ as outer, config, main
from antismash.common.errors import AntismashInputError
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
        comment = bio.annotations["structured_comment"]["antiSMASH-Data"]
        assert comment["Version"] == options.version
        assert "Original ID" not in comment
        assert "Starting at" not in comment
        assert "Ending at" not in comment

        bio.annotations["structured_comment"].pop("antiSMASH-Data")
        options.start = 7
        main.add_antismash_comments([(rec, bio)], options)
        comment = bio.annotations["structured_comment"]["antiSMASH-Data"]
        assert "Original ID" not in comment
        assert comment["Starting at"] == "7"

        bio.annotations["structured_comment"].pop("antiSMASH-Data")
        options.start = -1
        options.end = 1000
        main.add_antismash_comments([(rec, bio)], options)
        comment = bio.annotations["structured_comment"]["antiSMASH-Data"]
        assert "Original ID" not in comment
        assert comment["Ending at"] == "1000"

        bio.annotations["structured_comment"].pop("antiSMASH-Data")
        options.end = -1
        rec.original_id = "something else"
        main.add_antismash_comments([(rec, bio)], options)
        comment = bio.annotations["structured_comment"]["antiSMASH-Data"]
        assert comment["Original ID"] == rec.original_id

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


class TestErrorsInLogfile(unittest.TestCase):
    def setUp(self):
        destroy_config()

    def tearDown(self):
        destroy_config()

    # in a TestCase, the pytest capture system can't be used for a 'test_*' method
    # this solves the problem by storing it locally
    @pytest.fixture(autouse=True)
    def capture(self, capsys):
        self.outputs = capsys  # pylint: disable=attribute-defined-outside-init

    @patch.object(outer.os.path, "isfile", return_value=True)
    @patch.object(outer, "get_git_version", return_value="deadbeef")
    def test_error_logging_to_file(self, *_):
        error_message = "bad"
        with NamedTemporaryFile() as logfile:
            logfile_path = logfile.name
            # trick the input path checking by pointing it to the (empty) log file
            args = [logfile_path, "--logfile", logfile_path]
            options = build_config(args, isolated=True)
            with patch.object(main, "_run_antismash", side_effect=AntismashInputError(error_message)) as runner:
                # avoid constructing a non-isolated config
                with patch.object(config.AntismashParser, "parse_args", return_value=options):
                    assert build_config(args).logfile == logfile_path
                    assert get_config(args).logfile == logfile_path

                # finally, run the actual function being tested
                result = outer.main(args)

                # the run should count as a failure
                assert result == 1
                # the mocked inner call should have been run, with the options setting a logfile
                assert runner.call_count == 1
                assert runner.call_args == ((logfile_path, options),)
            # read the temporary log file before it's deleted
            with open(logfile_path, encoding="utf-8") as handle:
                log_lines = handle.read().splitlines()

        def check_message(line):
            assert line.startswith("ERROR") and line.endswith(error_message)

        # the error should be logged
        assert len(log_lines) == 1
        check_message(log_lines[0])

        # and also printed to stderr, since the logfile would hide it from the terminal
        printed = self.outputs.readouterr()  # fetches stdout/stderr via pytest
        check_message(printed.err.strip())
        # and nothing should have gone to stdout
        assert not printed.out


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
