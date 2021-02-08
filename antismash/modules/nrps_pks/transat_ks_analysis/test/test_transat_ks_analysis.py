# License: GNU Affero General Public License v3 or later
# A copy of GNU AGPL v3 should have been included in this software package in LICENSE.txt.

# for test files, silence irrelevant and noisy pylint warnings
# pylint: disable=no-self-use,protected-access,missing-docstring

import unittest

from antismash.modules.nrps_pks.data_structures import Prediction
from antismash.modules.nrps_pks.transat_ks_analysis.transat_ks_analysis import KSResult, KSPrediction, \
    get_leaf2clade, _LEAF2CLADE_TBL


class TestKSResult(unittest.TestCase):

    def setUp(self):
        self.result = KSResult("test_clade", "test_specificity", 0.0)

    def test_correct_instantiation(self):
        assert self.result.clade == "test_clade"
        assert self.result.specificity == "test_specificity"
        assert self.result.mass_score == 0.0

    def test_wrong_instantiation(self):
        self.assertRaises(AssertionError, KSResult, 1, "test_specificity", 0.0)
        self.assertRaises(AssertionError, KSResult, "test_clade", 1, 0.0)
        self.assertRaises(AssertionError, KSResult, "test_clade", "test_specificity", 0)

    def test_str(self):
        assert str(self.result) == "KSResult(clade=test_clade, specificity=test_specificity, mass_score=0.0)"

    def test_repr(self):
        assert repr(self.result) == "KSResult(clade=test_clade, specificity=test_specificity, mass_score=0.0)"

    def test_to_json(self):
        to_json_return = self.result.to_json()
        assert to_json_return == ("test_clade", "test_specificity", 0.0)

    def test_from_json(self):
        result = KSResult.from_json(("test_clade", "test_specificity", 0.0))
        assert result.clade == "test_clade"
        assert result.specificity == "test_specificity"
        assert result.mass_score == 0.0


class TestKSPrediction(unittest.TestCase):

    def setUp(self):
        results = {"test_specificity1": KSResult("test_clade1", "test_specificity1", 0.2),
                   "test_specificity2": KSResult("test_clade2", "test_specificity2", 0.7),
                   "test_specificity3": KSResult("test_clade3", "test_specificity3", 0.1)}
        self.prediction = KSPrediction(results)

    def test_correct_instantiation(self):
        # assert the order of predictions is as expected
        assert isinstance(self.prediction, Prediction)
        assert self.prediction.method == "transPACT"

        assert len(self.prediction.predictions) == 3
        assert self.prediction.predictions[0][0] == "test_specificity2"
        assert self.prediction.predictions[0][1].clade == "test_clade2"
        assert self.prediction.predictions[0][1].specificity == "test_specificity2"
        assert self.prediction.predictions[0][1].mass_score == 0.7

        assert self.prediction.predictions[1][0] == "test_specificity1"
        assert self.prediction.predictions[1][1].clade == "test_clade1"
        assert self.prediction.predictions[1][1].specificity == "test_specificity1"
        assert self.prediction.predictions[1][1].mass_score == 0.2

        assert self.prediction.predictions[2][0] == "test_specificity3"
        assert self.prediction.predictions[2][1].clade == "test_clade3"
        assert self.prediction.predictions[2][1].specificity == "test_specificity3"
        assert self.prediction.predictions[2][1].mass_score == 0.1

    def test_get_classification(self):
        empty_prediction = KSPrediction({})
        assert empty_prediction.get_classification() == []
        assert self.prediction.get_classification() == ["test_specificity2", "test_specificity1", "test_specificity3"]

    def test_as_html(self):
        empty_prediction = KSPrediction({})
        assert empty_prediction.as_html() == "No matches"
        assert self.prediction.as_html() == "<dl>\n " \
                                            "<dt>transPACT assigned specificiy:</dt>\n" \
                                            "<dd>test_specificity2 (test_clade2): 0.7 mass score (maximum=1.0)</dd>\n" \
                                            "<dd>test_specificity1 (test_clade1): 0.2 mass score (maximum=1.0)</dd>\n" \
                                            "<dd>test_specificity3 (test_clade3): 0.1 mass score (maximum=1.0)</dd>\n" \
                                            "</dl>\n"

    def test_to_json(self):
        to_json_return = self.prediction.to_json()
        assert to_json_return == {"method": "transPACT",
                                  "predictions": {"test_specificity2": ("test_clade2", "test_specificity2", 0.7),
                                                  "test_specificity1": ("test_clade1", "test_specificity1", 0.2),
                                                  "test_specificity3": ("test_clade3", "test_specificity3", 0.1)}}

    def test_from_json(self):
        prediction = KSPrediction.from_json(
            {"method": "transPACT", "predictions": {"test_specificity1": ("test_clade1", "test_specificity1", 0.2)}})
        assert prediction.method == "transPACT"
        assert len(prediction.predictions) == 1
        assert prediction.predictions[0][0] == "test_specificity1"
        assert prediction.predictions[0][1].clade == "test_clade1"
        assert prediction.predictions[0][1].specificity == "test_specificity1"
        assert prediction.predictions[0][1].mass_score == 0.2


class TestGetLeaf2Clade(unittest.TestCase):

    def test_leaf2clade_output(self):
        # test if the length of the output is as expected, more extensive testing seems unnecessary
        fun_clades, clade2ann = get_leaf2clade(_LEAF2CLADE_TBL)
        assert len(fun_clades) == 650
        assert len(clade2ann) == 104
