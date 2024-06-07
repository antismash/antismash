# License: GNU Affero General Public License v3 or later
# A copy of GNU AGPL v3 should have been included in this software package in LICENSE.txt.

# for test files, silence irrelevant and noisy pylint warnings
# pylint: disable=use-implicit-booleaness-not-comparison,protected-access,missing-docstring,consider-using-with

from argparse import Namespace
import glob
import os
from tempfile import NamedTemporaryFile, TemporaryDirectory
import unittest

import antismash
from antismash.main import run_antismash, get_all_modules, prepare_module_data
from antismash.common import path
from antismash.common.test.helpers import get_path_to_nisin_genbank
from antismash.common.secmet.test.helpers import rotate
from antismash.config import build_config, destroy_config, get_config, update_config


class TestAntismash(unittest.TestCase):
    def setUp(self):
        self.temp_dir = TemporaryDirectory()
        self.config = build_config(self.get_args() + ["--output-dir", self.temp_dir.name],
                                   isolated=True, modules=get_all_modules())
        self.input_file = get_path_to_nisin_genbank()
        self.file_prefix = "nisin"

    def get_args(self):
        return ["--minimal"]

    def tearDown(self):
        destroy_config()
        self.temp_dir.cleanup()

    def test_nisin_minimal(self):
        run_antismash(self.input_file, self.config)
        self.check_output_files()

    def check_output_files(self, filenames=None):
        out_dir = self.config.output_dir
        assert os.path.exists(out_dir)
        for filename in [f"{self.file_prefix}.json"] + (filenames or []):
            assert os.path.exists(os.path.join(out_dir, filename))


class TestSkipZip(TestAntismash):
    def get_args(self):
        return ["--minimal", "--skip-zip-file"]

    def check_output_files(self, filenames=None):
        super().check_output_files(filenames)
        assert not os.path.exists(os.path.join(self.config.output_dir, f"{self.file_prefix}.zip"))


class TestEnableHTML(TestAntismash):
    def get_args(self):
        return ["--minimal", "--enable-html"]

    def check_output_files(self, filenames=None):
        super().check_output_files(["index.html"])


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

    def check_output_files(self, filenames=None):
        super().check_output_files(filenames=[self.logfile.name])


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
        run_antismash(self.input_file, self.config)
        self.check_output_files()

        # remove html file and make sure it's recreated
        os.unlink(os.path.join(self.config.output_dir, "index.html"))
        update_config({"reuse_results": os.path.join(self.config.output_dir, f"{self.file_prefix}.json")})
        run_antismash(None, self.config)
        self.check_output_files()


class TestCircularReuse(TestResultsReuse):
    def setUp(self):
        super().setUp()
        record = antismash.common.secmet.Record.from_genbank(self.input_file)[0]
        rotate(record, len(record) // 2, padding=20_000)
        self.temp_file = NamedTemporaryFile(suffix=".gbk")
        self.input_file = self.temp_file.name
        self.file_prefix, _ = os.path.splitext(os.path.basename(self.input_file))
        record.to_genbank(self.input_file)

    def check_output_files(self, filenames=None):
        # start with typical checks
        super().check_output_files()

        out_gbk = os.path.join(self.config.output_dir, f"{self.file_prefix}.gbk")
        record = antismash.common.secmet.Record.from_genbank(out_gbk, taxon="bacteria")[0]
        assert record.is_circular()
        # and contains a protocluster over the origin
        protoclusters = record.get_protoclusters()
        assert len(protoclusters) == 1
        assert protoclusters[0].crosses_origin()


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

    def test_prereqs_with_container_data(self):
        modules = get_all_modules()
        options = build_config(["--databases", "/some/mounted_at_runtime/path"],
                               isolated=True, modules=modules)
        # there should be no errors for missing databases
        antismash.main.check_prerequisites(modules, options)

    def test_prepare_data_with_container_data(self):
        modules = get_all_modules()
        build_config(["--databases", "/some/mounted_at_runtime/path"],
                     isolated=True, modules=modules)
        # there should be no errors for missing databases
        prepare_module_data()

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
