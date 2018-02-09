# License: GNU Affero General Public License v3 or later
# A copy of GNU AGPL v3 should have been included in this software package in LICENSE.txt.

# for test files, silence irrelevant and noisy pylint warnings
# pylint: disable=no-self-use,protected-access,missing-docstring

import unittest

from antismash.common import utils


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
        self.assertEqual(802.9621, rpa.molecular_weight())  # default is True

        rpa = utils.RobustProteinAnalysis("MAGICXHAT", ignore_invalid=True)
        self.assertEqual(802.9621, rpa.molecular_weight())

    def test_molecular_weight_average(self):
        """ Test RobustProteinAnalysis.molecular_weight() calculates
            correct weight when not ignoring invalids
        """
        rpa = utils.RobustProteinAnalysis("MAGICXHAT", ignore_invalid=False)
        self.assertEqual(912.9621, rpa.molecular_weight())


class TestSignatureBuilding(unittest.TestCase):
    def test_extract_by_reference_positions(self):
        sig = utils.extract_by_reference_positions("ABC-DE-F", "A-BC-DEF", [0, 1, 3, 4])
        assert sig == "ACE-"
        sig = utils.extract_by_reference_positions("ABCDF", "ABCDE", [0, 1, 3, 4])
        assert sig == "ABDF"
