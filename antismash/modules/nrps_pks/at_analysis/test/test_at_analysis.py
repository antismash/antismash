# License: GNU Affero General Public License v3 or later
# A copy of GNU AGPL v3 should have been included in this software package in LICENSE.txt.

# for test files, silence irrelevant and noisy pylint warnings
# pylint: disable=use-implicit-booleaness-not-comparison,protected-access,missing-docstring

import unittest

from antismash.modules.nrps_pks.at_analysis import at_analysis


class TestScoring(unittest.TestCase):
    def test_score_threshold(self):
        ref_sigs = {f"{i}_{mon}": i*24 for i, mon in zip("ABC", ["mal", "mmal", "emal"])}
        query_sig = ["A"] * 24
        query_sig[3:7] = "B"
        query_sig[7:13] = "C"
        query_sigs = {"Q1": "".join(query_sig)}
        result = at_analysis.score_signatures(query_sigs, ref_sigs)
        assert len(result) == 1
        assert list(result) == ["Q1"]
        assert len(result["Q1"].predictions) == 1
        monomer, prediction = result["Q1"].predictions[0]
        assert monomer == "mal"
        assert prediction.name == "A_mal"
        assert prediction.signature == "A"*24
        self.assertAlmostEqual(prediction.score, 100*14/24)

    def test_hit_limit(self):
        ref_sigs = {str(i): "A"*24 for i in range(11)}
        query_sig = ["A"] * 24
        query_sigs = {"Q1": "".join(query_sig)}
        result = at_analysis.score_signatures(query_sigs, ref_sigs)
        assert len(result) == 1
        assert len(result["Q1"].predictions) == 11
        assert set(res.name for _, res in result["Q1"].predictions).issubset(set(ref_sigs))


class TestATPositions(unittest.TestCase):
    def test_structure(self):
        positions = at_analysis.get_at_positions()
        assert len(positions) == at_analysis._SIGNATURE_LENGTH
        for position in positions:
            assert isinstance(position, int)

    def test_start_position(self):
        positions = at_analysis.get_at_positions(startpos=0)
        assert positions[0] == 11
        positions = at_analysis.get_at_positions(startpos=7)
        assert positions[0] == 4
        positions = at_analysis.get_at_positions(startpos=9)
        assert positions[0] == 2
        positions = at_analysis.get_at_positions()
        assert positions[0] == 4


class TestATResult(unittest.TestCase):
    def test_bad_inputs(self):
        with self.assertRaises(AssertionError):
            at_analysis.ATResult("name", False, 0.5)
        with self.assertRaises(AssertionError):
            at_analysis.ATResult("name", "good", 0)
        with self.assertRaises(AssertionError):
            at_analysis.ATResult("name", "good", "str")
        with self.assertRaises(AssertionError):
            at_analysis.ATResult(False, "good", 0.5)

    def test_construction(self):
        result = at_analysis.ATResult("a", "b", 0.9)
        assert result.score == 0.9
        assert result.name == "a"
        assert result.signature == "b"

    def test_json(self):
        result = at_analysis.ATResult("a", "b", 0.9)
        regenned = at_analysis.ATResult.from_json(result.to_json())
        assert result.score == regenned.score
        assert result.name == regenned.name
        assert result.signature == regenned.signature
