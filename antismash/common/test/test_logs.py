# License: GNU Affero General Public License v3 or later
# A copy of GNU AGPL v3 should have been included in this software package in LICENSE.txt.

# for test files, silence irrelevant and noisy pylint warnings
# pylint: disable=no-self-use,protected-access,missing-docstring

import logging
import os
from tempfile import NamedTemporaryFile, TemporaryDirectory
import unittest

import pytest

from antismash.common.logs import changed_logging


class LogTest(unittest.TestCase):
    def setUp(self):
        self.logger = logging.getLogger()
        assert self.logger.getEffectiveLevel() == logging.WARNING

    # grab the pytest logging capture fixture
    @pytest.fixture(autouse=True)
    def inject_fixtures(self, caplog):
        caplog.set_level(logging.WARNING)
        self.caplog = caplog  # pylint: disable=attribute-defined-outside-init

    def test_debug(self):
        logging.debug("before")
        with changed_logging(debug=True):
            assert self.logger.getEffectiveLevel() == logging.DEBUG
            logging.debug("during")
        assert self.logger.getEffectiveLevel() == logging.WARNING
        logging.debug("after")
        for _, level, msg in self.caplog.record_tuples:
            assert level == logging.DEBUG
            assert "during" in msg

    def test_verbose(self):
        logging.info("before")
        with changed_logging(verbose=True):
            logging.debug("during")
            logging.info("during")
            assert self.logger.getEffectiveLevel() == logging.INFO
        logging.info("after")
        assert self.logger.getEffectiveLevel() == logging.WARNING
        for _, level, msg in self.caplog.record_tuples:
            assert level == logging.INFO
            assert "during" in msg

    def test_default(self):
        logging.warning("before")
        with changed_logging():
            logging.warning("during")
            logging.info("during")
            logging.debug("during")
            assert self.logger.getEffectiveLevel() == logging.WARNING
        logging.warning("after")
        assert self.logger.getEffectiveLevel() == logging.WARNING
        assert len(self.caplog.records) == 3
        for _, level, _ in self.caplog.record_tuples:
            assert level == logging.WARNING

    def test_debug_file(self):
        with NamedTemporaryFile() as logfile:
            with changed_logging(logfile=logfile.name, debug=True):
                logging.debug("during")
            logging.debug("after")

            with open(logfile.name) as handle:
                content = handle.read()
            assert "during" in content
            assert "after" not in content

    def test_logfile_always_verbose(self):
        with NamedTemporaryFile() as logfile:
            with changed_logging(logfile=logfile.name):
                logging.info("during")
                logging.debug("debug")
            logging.info("after")

            with open(logfile.name) as handle:
                content = handle.read()
            assert "during" in content
            assert "debug" not in content
            assert "after" not in content

    def test_logfile_default(self):
        with NamedTemporaryFile() as logfile:
            with changed_logging(logfile=logfile.name):
                logging.warning("during")
                logging.debug("debug")
            logging.warning("after")

            with open(logfile.name) as handle:
                content = handle.read()
            assert "during" in content
            assert "debug" not in content
            assert "after" not in content

    def test_logfile_debug(self):
        with NamedTemporaryFile() as logfile:
            with changed_logging(logfile=logfile.name, debug=True):
                logging.debug("during")
            logging.debug("after")

            with open(logfile.name) as handle:
                content = handle.read()
            assert "during" in content
            assert "after" not in content

    def test_logfile_nested_dirs(self):
        with TemporaryDirectory() as temp_dir:
            child = os.path.join(temp_dir, "child")
            assert not os.path.exists(child)
            logfile = os.path.join(child, "nested", "logfile")
            with changed_logging(logfile=logfile):
                logging.warning("during")
            assert os.path.exists(logfile)
