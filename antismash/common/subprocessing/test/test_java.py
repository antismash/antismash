# License: GNU Affero General Public License v3 or later
# A copy of GNU AGPL v3 should have been included in this software package in LICENSE.txt.

# for test files, silence irrelevant and noisy pylint warnings
# pylint: disable=no-self-use,protected-access,missing-docstring

import unittest
from unittest.mock import patch

import antismash
from antismash.config import build_config, destroy_config

from antismash.common import subprocessing

class DummyResult(subprocessing.RunResult):
    def __init__(self, stderr: str):
        super().__init__(["dummy"], b"", stderr.encode(), 0, True, True)


OPENJDK_VERSION = """openjdk version "11.0.6" 2020-01-14
OpenJDK Runtime Environment (build 11.0.6+10-post-Ubuntu-1ubuntu118.04.1)
OpenJDK 64-Bit Server VM (build 11.0.6+10-post-Ubuntu-1ubuntu118.04.1, mixed mode, sharing)
"""

JAVA_VERSION = """java version "9.0.4"
Java(TM) SE Runtime Environment (build 9.0.4+11)
Java HotSpot(TM) 64-Bit Server VM (build 9.0.4+11, mixed mode)
"""

WITH_JAVA_OPTIONS = """Picked up _JAVA_OPTIONS: -Dawt.useSystemAAFontSettings=on
openjdk version "11.0.6" 2020-01-14
OpenJDK Runtime Environment (build 11.0.6+10-post-Ubuntu-1ubuntu118.04.1)
OpenJDK 64-Bit Server VM (build 11.0.6+10-post-Ubuntu-1ubuntu118.04.1, mixed mode, sharing)
"""

INVALID = """Hi there.
I'm some invalid "java" output.
I should trigger an error.
"""

class TestJava(unittest.TestCase):

    def setUp(self):
        build_config([], isolated=True, modules=antismash.get_all_modules())

    def tearDown(self):
        destroy_config()

    @patch("antismash.common.subprocessing.java.execute", return_value=DummyResult(OPENJDK_VERSION))
    def test_java_version_openjdk(self, _):
        version = subprocessing.run_java_version()
        assert version == "11.0.6"

    @patch("antismash.common.subprocessing.java.execute", return_value=DummyResult(WITH_JAVA_OPTIONS))
    def test_java_version_openjdk_options(self, _):
        version = subprocessing.run_java_version()
        assert version == "11.0.6"

    @patch("antismash.common.subprocessing.java.execute", return_value=DummyResult(JAVA_VERSION))
    def test_java_version_java(self, _):
        version = subprocessing.run_java_version()
        assert version == "9.0.4"

    @patch("antismash.common.subprocessing.java.execute", return_value=DummyResult(INVALID))
    def test_java_version_invalid(self, _):
        with self.assertRaises(RuntimeError):
            subprocessing.run_java_version()
