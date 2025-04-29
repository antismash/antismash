# License: GNU Affero General Public License v3 or later
# A copy of GNU AGPL v3 should have been included in this software package in LICENSE.txt.

# for test files, silence irrelevant and noisy pylint warnings
# pylint: disable=use-implicit-booleaness-not-comparison,protected-access,missing-docstring

import os
import time
import unittest
from unittest.mock import patch

from helperlibs.wrappers.io import TemporaryDirectory

from antismash.common import subprocessing, path
from antismash.config import build_config, destroy_config, get_config, update_config
from antismash.main import get_all_modules


def dummy(value=None):
    """ Since functions passed to parallel_function can't be locally defined,
        this has to be declared here
    """
    time.sleep(.1)
    if value == 1:
        raise ValueError("Lucky number")
    return os.getpid()


def busy_wait():
    """ Another function to pass to parallel_function, this time with a busy
        wait to avoid potentially parking multiple processes on the same core
        due to sleeping
    """
    start = time.time()
    while (time.time() - start) < 1:
        # add a little busy work so that time isn't called so often
        assert "X" not in "A" * 100000
    return os.getpid()


class TestParallelPython(unittest.TestCase):
    def setUp(self):
        self.config_cpus = 2
        build_config(["--cpus", str(self.config_cpus)])
        assert get_config().executables
        assert get_config().cpus == 2

    def tearDown(self):
        destroy_config()

    def test_dummy(self):
        # checks that the dummy function actually works as expected
        assert dummy() == os.getpid()
        assert dummy(0) == os.getpid()
        with self.assertRaisesRegex(ValueError, "Lucky number"):
            dummy(1)
        for i in range(2, 10):
            assert dummy(i) == os.getpid()

    def test_busy_wait(self):
        assert busy_wait() == os.getpid()

    def test_actually_parallel(self):
        results = subprocessing.parallel_function(busy_wait, [[]]*2, timeout=5)
        assert len(set(results)) == self.config_cpus
        # since it uses a pool, if we try with more we'll reuse those some procs
        results = subprocessing.parallel_function(busy_wait, [[]]*4, timeout=10)
        assert len(set(results)) == self.config_cpus
        # but check it did run the requested amount of times
        assert len(results) == 4

    def test_cpu_param(self):
        results = subprocessing.parallel_function(dummy, [[]]*4, cpus=4)
        assert len(set(results)) == 4
        assert len(results) == 4

    def test_timeout(self):
        start = time.time()
        timeout = 1
        with self.assertRaisesRegex(RuntimeError, "Timeout in parallel .*"):
            subprocessing.parallel_function(time.sleep, [[5]]*2, timeout=timeout)
        elapsed = time.time() - start
        assert elapsed < 1.5

    def test_results(self):
        # test single-arg functions
        res = subprocessing.parallel_function(len, [["a"], ["aa"], ["aaa"]])
        assert res == [1, 2, 3]

        # test multi-arg functions
        res = subprocessing.parallel_function(zip, [["%d" % i, "b"] for i in range(3)])
        assert [list(i) for i in res] == [[(str(i), "b")] for i in range(3)]

    def test_errors_propagated(self):
        with self.assertRaisesRegex(ValueError, "Lucky number"):
            subprocessing.parallel_function(dummy, [[i] for i in range(5)])

    def test_local_funcs(self):
        # this test is mostly to let us know if something changes in pool/pickle
        def local(value):
            return value + 1
        # ensure the function works as expected when called directly
        assert local(1) == 2
        # check if it still fails within a parallel pool
        with self.assertRaisesRegex(AttributeError, r"Can't (\w+) local object"):
            subprocessing.parallel_function(local, [[i] for i in range(3)])


class TestExecute(unittest.TestCase):
    def setUp(self):
        build_config([])

    def tearDown(self):
        destroy_config()

    def test_piping(self):
        result = subprocessing.execute(["pwd"])
        assert result.stdout.strip() == os.getcwd()
        assert result.stderr.strip() == ""
        assert not result.return_code and result.successful()

        result = subprocessing.execute(["cat", "--bad-option"])
        assert result.stdout.strip() == ""
        # the exact message isn't important, only that cat emitted it
        assert result.stderr.startswith("cat: ")
        assert result.return_code and not result.successful()

        result = subprocessing.execute(["cat"], stdin="fish")
        assert result.stdout.strip() == "fish"
        assert not result.stderr
        assert not result.return_code and result.successful()

    def test_redirection(self):
        result = subprocessing.execute(["echo", "test"], stdout=open(os.devnull, "w"))
        with self.assertRaisesRegex(ValueError, "stdout was redirected to file, unable to access"):
            result.stdout.strip()
        assert not result.stderr
        assert not result.return_code
        assert result.get_command_string() == "echo test"

        result = subprocessing.execute(["cat", "--bad-option"], stderr=open(os.devnull, "w"))
        assert result.stdout.strip() == ""
        with self.assertRaisesRegex(ValueError, "stderr was redirected to file, unable to access"):
            _ = result.stderr
        assert result.return_code

    def test_timeout(self):
        start = time.time()
        with self.assertRaisesRegex(RuntimeError, "Child process .* timed out after 1 second"):
            subprocessing.execute(["sleep", "10"], timeout=1)
        elapsed = time.time() - start
        assert elapsed < 1.5

    def test_env(self):
        class Dummy:
            def __init__(self):
                self.returncode = 0

            def communicate(self, *_args, **_kwargs):
                return ("stderr".encode(), "stdout".encode())

            def __enter__(self, *_args, **_kwargs):
                return Dummy()

            def __exit__(self, *_args, **_kwargs):
                pass

        dummy_proc = Dummy()
        with patch.object(subprocessing.base, "Popen", return_value=dummy_proc) as patched:
            subprocessing.execute(["dummy"], environment_overrides={"something": "value"})
            assert patched.call_args.kwargs["env"]["something"] == "value"

            patched.reset()

            subprocessing.execute(["dummy"], environment_overrides={"PATH": "/some/dummy/path"})
            assert patched.call_args.kwargs["env"]["PATH"] == "/some/dummy/path"


class TestParallelExecute(unittest.TestCase):
    def setUp(self):
        self.config_cpus = 2
        build_config(["--cpus", str(self.config_cpus)])
        assert get_config().cpus == 2

    def tearDown(self):
        destroy_config()

    def test_actually_parallel(self):
        start = time.time()
        subprocessing.parallel_execute([["sleep", "1"]]*2)
        elapsed = time.time() - start
        assert elapsed < 1.5

    def test_timeout(self):
        start = time.time()
        with self.assertRaisesRegex(RuntimeError, ""):
            subprocessing.parallel_execute([["sleep", "10"]]*2, timeout=1)
        elapsed = time.time() - start
        assert elapsed < 1.5


class TestDiamondDatabaseChecks(unittest.TestCase):
    def setUp(self):
        self.format0_file = path.get_full_path(__file__, "data", "format0.dmnd")
        self.format1_file = path.get_full_path(__file__, "data", "format1.dmnd")
        self.empty = path.get_full_path(__file__, "data", "empty.dmnd")

        options = build_config([], isolated=True, modules=get_all_modules())
        self.old_config = get_config().__dict__
        self.options = update_config(options)

    def tearDown(self):
        destroy_config()
        update_config(self.old_config)

    def test_check_diamond_db_compatible(self):
        with TemporaryDirectory(change=True):
            dummy_fasta = "dummy.fa"
            dummy_db = "dummy.dmnd"
            with open(dummy_fasta, "w") as handle:
                handle.write(">test\nM\n")
            subprocessing.run_diamond_makedb(dummy_db, dummy_fasta)
            compatible_format = subprocessing.diamond._extract_db_format(dummy_db)
            assert subprocessing.diamond.check_diamond_db_compatible(dummy_db)

        broken_file = self.format0_file if compatible_format > 0 else self.format1_file

        assert not subprocessing.diamond.check_diamond_db_compatible(broken_file)
        assert not subprocessing.diamond.check_diamond_db_compatible(self.empty)
