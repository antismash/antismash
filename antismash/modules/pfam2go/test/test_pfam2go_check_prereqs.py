# License: GNU Affero General Public License v3 or later
# A copy of GNU AGPL v3 should have been included in this software package in LICENSE.txt.

# for test files, silence irrelevant and noisy pylint warnings
# pylint: disable=use-implicit-booleaness-not-comparison,protected-access,missing-docstring

import unittest
from unittest.mock import patch

from antismash.common import path
from antismash.modules.pfam2go import check_prereqs


class Pfam2GoPrereqsTest(unittest.TestCase):
    def test_check_prereqs(self):
        ret = check_prereqs(None)
        assert ret == []

    def test_check_missing_file(self):
        with patch.object(path, 'locate_file', return_value=None):
            ret = check_prereqs(None)
        assert 'Failed to locate Pfam to Gene Ontology mapping file' in ret
