# License: GNU Affero General Public License v3 or later
# A copy of GNU AGPL v3 should have been included in this software package in LICENSE.txt.

# for test files, silence irrelevant and noisy pylint warnings
# pylint: disable=use-implicit-booleaness-not-comparison,protected-access,missing-docstring

import unittest

from antismash.config import destroy_config, build_config
from antismash.modules.lassopeptides import check_prereqs, local_config


class TestCheckPrereqs(unittest.TestCase):
    def setUp(self):
        self.options = build_config([])

    def tearDown(self):
        destroy_config()

    def test_check_prereqs(self):
        "Test lassopeptides.check_prereqs()"
        ret = check_prereqs(self.options)
        self.assertEqual(ret, [])
        self.assertTrue(local_config().fimo_present)
