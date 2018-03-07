# License: GNU Affero General Public License v3 or later
# A copy of GNU AGPL v3 should have been included in this software package in LICENSE.txt.

# for test files, silence irrelevant and noisy pylint warnings
# pylint: disable=no-self-use,protected-access,missing-docstring

import os
import tempfile
import unittest

from antismash.common import pfamdb


class TestPFAMs(unittest.TestCase):
    def setUp(self):
        pfamdb.KNOWN_MAPPINGS.clear()

    def test_version_finder(self):
        for path, version in [("/path/to/pfam/31.0/bits.hmm", "31.0"),
                              ("relpath/to/pfam/31.0", "31.0"),
                              ("/path/to/other/pfam/5.1/stuff", "5.1")]:
            assert pfamdb.get_db_version_from_path(path) == version
        for path in ["/bad/path/31.0", "/other/bad/pfams/31.0"]:
            with self.assertRaisesRegex(ValueError, "Database path does not contain a 'pfam' directory"):
                pfamdb.get_db_version_from_path(path)
        for path in ["/bad/path/pfam/31.0x", "/other/bad/pfam/31.0x/stuff"]:
            with self.assertRaisesRegex(ValueError, "No valid database version found in PFAM database path"):
                pfamdb.get_db_version_from_path(path)

    def test_find_latest(self):
        with tempfile.TemporaryDirectory(prefix="aS.pfamdbtest") as temp_db_layout:
            os.makedirs(os.path.join(temp_db_layout, "pfam", "30.7invalid"))
            os.makedirs(os.path.join(temp_db_layout, "pfam", "invalid30.7"))
            os.makedirs(os.path.join(temp_db_layout, "pfam", "irrelevant"))

            with self.assertRaisesRegex(Exception, "No matching PFAM database in location " + temp_db_layout):
                pfamdb.find_latest_database_version(temp_db_layout)

            os.makedirs(os.path.join(temp_db_layout, "pfam", "31.0"))
            assert pfamdb.find_latest_database_version(temp_db_layout) == "31.0"
            os.makedirs(os.path.join(temp_db_layout, "pfam", "31.2"))
            assert pfamdb.find_latest_database_version(temp_db_layout) == "31.2"
            os.makedirs(os.path.join(temp_db_layout, "pfam", "30.7"))
            assert pfamdb.find_latest_database_version(temp_db_layout) == "31.2"
