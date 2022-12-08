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
            ("method_a", "Ala")
        ]
        assert "Ala" == self.go(_generate_predictions(data))

    def test_majority(self):
        data = [
            ("method_a", "Ala"),
            ("method_b", "Gly"),
            ("method_c", "Ala"),
        ]
        assert "Ala" == self.go(_generate_predictions(data))

    def test_tie(self):
        data = [
            ("method_a", "Ala"),
            ("method_b", "Gly"),
        ]
        assert "X" == self.go(_generate_predictions(data))

    def test_tie_with_extra_predictions(self):
        data = [
            ("method_a", "Ala"),
            ("method_b", "Gly"),
            ("method_c", "Ala"),
            ("method_d", "Gly"),
            ("method_e", "Tyr"),
        ]
        assert "X" == self.go(_generate_predictions(data))

def _generate_predictions(data: list[tuple[str, str]]) -> dict[str, Prediction]:
    return {
        method: SimplePrediction(method, prediction) for method, prediction in data
    }
