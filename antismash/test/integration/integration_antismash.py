# License: GNU Affero General Public License v3 or later
# A copy of GNU AGPL v3 should have been included in this software package in LICENSE.txt.

# for test files, silence irrelevant and noisy pylint warnings
# pylint: disable=no-self-use,protected-access,missing-docstring,consider-using-with

from argparse import Namespace
import glob
import os
from tempfile import NamedTemporaryFile, TemporaryDirectory
import unittest

import antismash
from antismash.main import run_antismash, get_all_modules, prepare_module_data
from antismash.common import path
from antismash.common.test.helpers import get_path_to_nisin_genbank
from antismash.config import build_config, destroy_config, get_config, update_config


class TestAntismash(unittest.TestCase):
    def setUp(self):
        self.temp_dir = TemporaryDirectory()
        self.config = build_config(self.get_args() + ["--output-dir", self.temp_dir.name],
                                   isolated=True, modules=get_all_modules())

    def get_args(self):
        return ["--minimal"]

    def tearDown(self):
        destroy_config()
        self.temp_dir.cleanup()

    def test_nisin_minimal(self):
        run_antismash(get_path_to_nisin_genbank(), self.config)
        self.check_output_files()

    def check_output_files(self):
        out_dir = self.config.output_dir
        assert os.path.exists(out_dir)
        for filename in ["nisin.json"]:
            assert os.path.exists(os.path.join(out_dir, filename))


class TestSkipZip(TestAntismash):
    def get_args(self):
        return ["--minimal", "--skip-zip-file"]

    def check_output_files(self):
        super().check_output_files()
        assert not os.path.exists(os.path.join(self.config.output_dir, "nisin.zip"))


class TestEnableHTML(TestAntismash):
    def get_args(self):
        return ["--minimal", "--enable-html"]

    def check_output_files(self):
        super().check_output_files()
        assert os.path.exists(os.path.join(self.config.output_dir, "index.html"))


class TestProfiling(TestAntismash):
    def get_args(self):
        return ["--minimal", "--profiling"]


class TestDebug(TestAntismash):
    def get_args(self):
        return ["--minimal", "-d"]


class TestVerbose(TestAntismash):
    def get_args(self):
        return ["--minimal", "-v"]


class TestLogging(TestAntismash):
    def setUp(self):
        self.logfile = NamedTemporaryFile()
        super().setUp()

    def tearDown(self):
        super().tearDown()
        self.logfile.close()

    def get_args(self):
        return ["--minimal", "-v", "--logfile", os.path.join(self.logfile.name)]

    def check_output_files(self):
        super().check_output_files()
        assert os.path.exists(self.logfile.name)


class TestResultsReuse(TestAntismash):
    def get_args(self):
        return ["--minimal", "--enable-html"]

    def test_nisin_minimal(self):
        # make sure the output directory isn't filled
        out_dir = self.config.output_dir
        assert not list(glob.glob(os.path.join(out_dir, "*")))

        # die with neither inputs provided
        with self.assertRaisesRegex(ValueError, "No sequence file or prior results to read"):
            run_antismash(None, self.config)

        # make sure no files created
        assert not list(glob.glob(os.path.join(out_dir, "*")))

        # do a normal run
        run_antismash(get_path_to_nisin_genbank(), self.config)
        self.check_output_files()

        # remove html file and make sure it's recreated
        os.unlink(os.path.join(self.config.output_dir, "index.html"))
        update_config({"reuse_results": os.path.join(self.config.output_dir, "nisin.json")})
        run_antismash(None, self.config)
        self.check_output_files()


class TestModuleData(unittest.TestCase):
    def setUp(self):
        build_config([], isolated=True, modules=get_all_modules())

    def tearDown(self):
        destroy_config()

    def test_check_prereqs(self):
        options = build_config(["--check-prereqs"], isolated=False, modules=get_all_modules())
        assert run_antismash("", options) == 0

    def test_check_prereqs_missing_executables(self):
        options = build_config(["--check-prereqs"], isolated=True, modules=get_all_modules())
        update_config({"executables": Namespace()})
        assert hasattr(get_config(), "executables")
        assert not get_config().executables.__dict__
        with self.assertRaisesRegex(RuntimeError, "failing prereq"):
            antismash.main.check_prerequisites(get_all_modules(), options)

    def test_prepare_module_data(self):
        # make sure there's some to start with
        search = path.get_full_path(antismash.__file__, '**', "*.h3?")
        existing_press_files = glob.glob(search, recursive=True)
        assert existing_press_files

        # then remove them all
        for pressed in existing_press_files:
            os.unlink(pressed)
        current_press_files = glob.glob(search, recursive=True)
        assert not current_press_files

        # and make sure they're all regenerated properly
        prepare_module_data()
        current_press_files = glob.glob(search, recursive=True)
        assert current_press_files == existing_press_files
