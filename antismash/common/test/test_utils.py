# License: GNU Affero General Public License v3 or later
# A copy of GNU AGPL v3 should have been included in this software package in LICENSE.txt.

# for test files, silence irrelevant and noisy pylint warnings
# pylint: disable=no-self-use,protected-access,missing-docstring

import unittest

from antismash.common import utils

class TestUniqueID(unittest.TestCase):
    def test_bad_starts(self):
        for bad_start in ["start", None, dict(), list()]:
            with self.assertRaises((ValueError, TypeError)):
                utils.generate_unique_id("pref", [], bad_start)

    def test_bad_collections(self):
        class Dummy:
            pass
        for bad_existing in [None, Dummy(), 2]:
            with self.assertRaises((ValueError, TypeError)):
                utils.generate_unique_id("pref", bad_existing, 1)

    def test_bad_max(self):
        for bad_max in ["start", None, dict(), list()]:
            with self.assertRaises((ValueError, TypeError)):
                utils.generate_unique_id("pref", {}, 1, bad_max)

    def test_generation(self):
        existing = {"a_%d" % i for i in range(15)}
        new, counter = utils.generate_unique_id("a", existing)
        assert len(existing) == 15 and new not in existing
        assert new == "a_15" and counter == 15

        new, counter = utils.generate_unique_id("a", existing, start=17)
        assert len(existing) == 15 and new not in existing
        assert new == "a_17" and counter == 17

        new, counter = utils.generate_unique_id("b", existing)
        assert len(existing) == 15 and new not in existing
        assert new == "b_0" and counter == 0


    def test_overlong(self):
        existing = {"a_%d" % i for i in range(150)}
        # prefix itself too long
        with self.assertRaisesRegex(RuntimeError, "Could not generate .*"):
            utils.generate_unique_id("aaa", existing, start=0, max_length=3)

        # the generated number is too long
        with self.assertRaisesRegex(RuntimeError, "Could not generate .*"):
            utils.generate_unique_id("a", existing, start=140, max_length=4)

class TestRobustProteinAnalysis(unittest.TestCase):
    def test_init(self):
        """Test RobustProteinAnalysis initialisation"""
        rpa = utils.RobustProteinAnalysis("MAGICHAT", ignore_invalid=True)
        self.assertIsInstance(rpa, utils.RobustProteinAnalysis)

        rpa = utils.RobustProteinAnalysis("MAGICHAT", ignore_invalid=False)
        self.assertIsInstance(rpa, utils.RobustProteinAnalysis)

        for bad_invalid in ["none", None, 3, []]:
            with self.assertRaises(TypeError):
                utils.RobustProteinAnalysis("MAGICHAT", ignore_invalid=bad_invalid)

    def test_uppercase(self):
        """Test RobustProteinAnalysis converts passed sequence to upper case"""
        rpa = utils.RobustProteinAnalysis("Magichat")
        assert rpa.original_sequence == "MAGICHAT"
        assert rpa.sequence == "MAGICHAT"

    def test_molecular_weight_ignore(self):
        """ Test RobustProteinAnalysis.molecular_weight() calculates
            correct weight when ignoring invalids"""
        rpa = utils.RobustProteinAnalysis("MAGICXHAT")
        self.assertEqual(802.9621, rpa.molecular_weight()) # default is True

        rpa = utils.RobustProteinAnalysis("MAGICXHAT", ignore_invalid=True)
        self.assertEqual(802.9621, rpa.molecular_weight())

    def test_molecular_weight_average(self):
        """ Test RobustProteinAnalysis.molecular_weight() calculates
            correct weight when not ignoring invalids
        """
        rpa = utils.RobustProteinAnalysis("MAGICXHAT", ignore_invalid=False)
        self.assertEqual(912.9621, rpa.molecular_weight())
