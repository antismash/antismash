# License: GNU Affero General Public License v3 or later
# A copy of GNU AGPL v3 should have been included in this software package in LICENSE.txt.

import unittest
from minimock import Mock, mock, restore, TraceTracker, assert_same_trace
# pylint: disable=unused-import
import antismash # used in mocking
# pylint: enable=unused-import
from antismash.modules.lanthipeptides import check_prereqs, get_config

class TestCheckPrereqs(unittest.TestCase):
    def setUp(self):
        self.maxDiff = None
        self.tt = TraceTracker()
        self.locate_exe = Mock('antismash.common.path.locate_executable',
                               tracker=self.tt, returns="/fake/path/to/binary")
        mock('antismash.common.path.locate_executable',
             mock_obj=self.locate_exe, tracker=self.tt)

    def tearDown(self):
        restore()

    def test_check_prereqs(self):
        "Test lanthipeptides.check_prereqs()"
        ret = check_prereqs()
        self.assertEqual(ret, [])
        expected = """    Called antismash.common.path.locate_executable('hmmpfam2')
    Called antismash.common.path.locate_executable('fimo')"""
        assert_same_trace(self.tt, expected)
        self.assertTrue(get_config().fimo_present)

    def test_check_binary_prereqs_failing(self):
        "Test lanthipeptides.check_prereqs() returns 'missing binary' error"
        self.locate_exe.mock_returns = None
        ret = check_prereqs()
        self.assertEqual(len(ret), 1)
        self.assertIn("Failed to locate executable for 'hmmpfam2'", ret)
        self.assertFalse(get_config().fimo_present)
