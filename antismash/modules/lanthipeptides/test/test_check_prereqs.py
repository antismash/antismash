# License: GNU Affero General Public License v3 or later
# A copy of GNU AGPL v3 should have been included in this software package in LICENSE.txt.

# for test files, silence irrelevant and noisy pylint warnings
# pylint: disable=use-implicit-booleaness-not-comparison,protected-access,missing-docstring

import unittest

from antismash.config import build_config, destroy_config
from antismash.modules.lanthipeptides import check_prereqs, get_config


class TestCheckPrereqs(unittest.TestCase):
    def setUp(self):
        self.config = build_config([])

    def tearDown(self):
        destroy_config()

    def test_check_prereqs(self):
        "Test lanthipeptides.check_prereqs()"
        ret = check_prereqs(self.config)
        self.assertEqual(ret, [])
        self.assertTrue(get_config().fimo_present)

    def test_check_binary_prereqs_failing(self):
        "Test lanthipeptides.check_prereqs() returns 'missing binary' error"
        self.config.executables.__dict__.pop("hmmpfam2")
        self.config.executables.__dict__.pop("fimo")
        ret = check_prereqs(self.config)
        self.assertEqual(len(ret), 1)
        self.assertIn("Failed to locate executable for 'hmmpfam2'", ret)
        self.assertFalse(get_config().fimo_present)
