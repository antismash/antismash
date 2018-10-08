# License: GNU Affero General Public License v3 or later
# A copy of GNU AGPL v3 should have been included in this software package in LICENSE.txt.

# for test files, silence irrelevant and noisy pylint warnings
# pylint: disable=no-self-use,protected-access,missing-docstring

import glob
import os
from tempfile import NamedTemporaryFile, TemporaryDirectory
import unittest

import antismash
from antismash.main import run_antismash, get_all_modules, prepare_module_data
from antismash.common import path
from antismash.common.test.helpers import get_path_to_nisin_genbank
from antismash.config import build_config, update_config, destroy_config


class TestAntismash(unittest.TestCase):
    def setUp(self):
        self.temp_dir = TemporaryDirectory()
        self.build_config(self.get_args())

    def get_args(self):
        return ["--minimal"]

    def build_config(self, args):
        self.default_options = build_config(args, isolated=True, modules=get_all_modules())
        self.default_options.all_enabled_modules = []
        self.default_options.output_dir = self.temp_dir.name

    def tearDown(self):
        destroy_config()
        self.temp_dir.cleanup()

    def test_nisin_minimal(self):
        run_antismash(get_path_to_nisin_genbank(), self.default_options)
        self.check_output_files()

    def check_output_files(self):
        out_dir = self.default_options.output_dir
        assert os.path.exists(out_dir)
        for filename in ["nisin.json", "index.html"]:
            assert os.path.exists(os.path.join(out_dir, filename))


class TestSkipZip(TestAntismash):
    def get_args(self):
        return ["--minimal", "--skip-zip-file"]

    def check_output_files(self):
        super().check_output_files()
        assert not os.path.exists(os.path.join(self.default_options.output_dir, "nisin.zip"))


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
    def test_nisin_minimal(self):
        # make sure the output directory isn't filled
        out_dir = self.default_options.output_dir
        assert not list(glob.glob(os.path.join(out_dir, "*")))

        # die with neither inputs provided
        with self.assertRaisesRegex(ValueError, "No sequence file or prior results to read"):
            run_antismash(None, self.default_options)

        # make sure no files created
        assert not list(glob.glob(os.path.join(out_dir, "*")))

        # do a normal run
        run_antismash(get_path_to_nisin_genbank(), self.default_options)
        self.check_output_files()

        # remove html file and make sure it's recreated
        os.unlink(os.path.join(self.default_options.output_dir, "index.html"))
        update_config({"reuse_results": os.path.join(self.default_options.output_dir, "nisin.json")})
        run_antismash(None, self.default_options)
        self.check_output_files()


class TestModuleData(unittest.TestCase):
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
