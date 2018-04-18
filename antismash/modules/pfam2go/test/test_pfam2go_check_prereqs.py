# License: GNU Affero General Public License v3 or later
# A copy of GNU AGPL v3 should have been included in this software package in LICENSE.txt.

# for test files, silence irrelevant and noisy pylint warnings
# pylint: disable=no-self-use,protected-access,missing-docstring

import unittest

from minimock import mock, Mock, restore

from antismash.common import path  # mocked, pylint: disable=unused-import
from antismash.modules.pfam2go import check_prereqs


class Pfam2GoPrereqsTest(unittest.TestCase):
    def setUp(self):
        self.locate_file = Mock('antismash.common.path.locate_file', returns="/fake/path/to/file")
        mock('path.locate_file', mock_obj=self.locate_file)

    def tearDown(self):
        restore()

    def test_check_prereqs(self):
        ret = check_prereqs()
        assert ret == []

    def test_check_missing_file(self):
        self.locate_file.mock_returns = None
        ret = check_prereqs()
        assert 'Failed to locate Pfam to Gene Ontology mapping file' in ret
