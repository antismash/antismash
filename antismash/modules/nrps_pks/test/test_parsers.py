# License: GNU Affero General Public License v3 or later
# A copy of GNU AGPL v3 should have been included in this software package in LICENSE.txt.

# for test files, silence irrelevant and noisy pylint warnings
# pylint: disable=use-implicit-booleaness-not-comparison,protected-access,missing-docstring

import unittest

from antismash.modules.nrps_pks import parsers
from antismash.modules.nrps_pks.data_structures import Prediction, SimplePrediction

class TestNrpsConsensus(unittest.TestCase):
    def setUp(self) -> None:
        self.go = parsers.generate_nrps_consensus

    def test_single(self):
        data = [
            ("method_a", "Ala", None)
        ]
        assert "Ala" == self.go(_generate_predictions(data))

    def test_majority(self):
        data = [
            ("method_a", "Ala", None),
            ("method_b", "Gly", None),
            ("method_c", "Ala", None),
        ]
        assert "Ala" == self.go(_generate_predictions(data))

    def test_tie(self):
        data = [
            ("method_a", "Ala", None),
            ("method_b", "Gly", None),
        ]
        assert "X" == self.go(_generate_predictions(data))

    def test_tie_with_extra_predictions(self):
        data = [
            ("method_a", "Ala", None),
            ("method_b", "Gly", None),
            ("method_c", "Ala", None),
            ("method_d", "Gly", None),
            ("method_e", "Tyr", None),
        ]
        assert "X" == self.go(_generate_predictions(data))

    def test_norine_differs(self):
        data = [
            ("method_a", "hydrophobic-aliphatic", "X")
        ]
        assert "X" == self.go(_generate_predictions(data))

def _generate_predictions(data: list[tuple[str, str, str]]) -> dict[str, Prediction]:
    return {
        method: SimplePrediction(method, prediction, norine) for method, prediction, norine in data
    }
