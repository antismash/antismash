import unittest
from math import isclose

from markupsafe import Markup

from antismash.modules.nrps_pks.paras import ParasModel, GeneralParasResult, ParasResult, ParasectResult


class MockMinimalResult:
    def __init__(self, signature, prediction_labels, predictions):
        self.signature = signature
        self.prediction_labels = prediction_labels
        self.predictions = predictions


class TestGeneralParasResult(unittest.TestCase):
    def setUp(self):
        self.paras_result = GeneralParasResult(ParasModel.PARAS,
                                               "ILIKEDATAEVENFAKEDATADIDISAYILIKED",
                                               ["tryptophan", "tyrosine", "phenylalanine"],
                                               [0.9, 0.05, 0.05],
                                               2)
        self.parasect_result = GeneralParasResult(ParasModel.PARASECT,
                                                  "ILIKEDATAEVENFAKEDATADIDISAYILIKED",
                                                  ["tryptophan", "tyrosine", "phenylalanine"],
                                                  [0.9, 0.8, 0.6],
                                                  3)

    def test_get_classification(self):
        assert self.paras_result.get_classification() == ["Trp"]
        assert self.parasect_result.get_classification() == ["Trp", "Tyr"]

    def test_classification_uncertain(self):
        pred = GeneralParasResult(ParasModel.PARAS,
                                  "ILIKEDATAEVENFAKEDATADIDISAYILIKED",
                                  ["tryptophan", "tyrosine", "phenylalanine"],
                                  [0.1, 0.05, 0.05],
                                  2)
        assert pred.uncertain
        pred = GeneralParasResult(ParasModel.PARASECT,
                                  "ILIKEDATAEVENFAKEDATADIDISAYILIKED",
                                  ["tryptophan", "tyrosine", "phenylalanine"],
                                  [0.6, 0.5, 0.5],
                                  2)

        assert pred.uncertain

    def test_classification_certain(self):
        pred = self.paras_result
        assert not pred.uncertain

        pred = self.parasect_result
        assert not pred.uncertain

    def test_valid(self):

        pred = self.paras_result

        assert pred.aa34 == "ILIKEDATAEVENFAKEDATADIDISAYILIKED"
        assert pred.method == "paras"
        assert pred.predicted_substrate == "tryptophan"
        assert isclose(pred.top_confidences[0], 0.9)
        assert pred.top_substrates == ["tryptophan", "tyrosine"]
        self.assertNearlyEqual(pred.top_confidences, [0.9, 0.05, 0.05])
        assert not pred.uncertain

        pred = self.parasect_result
        assert pred.method == "parasect"
        assert not pred.uncertain

    def test_invalid(self):
        with self.assertRaises(AssertionError):
            GeneralParasResult(ParasModel.PARAS,
                               "ILIKEDATAEVENFAKEDATADIDISAYILIKED",
                               ["tryptophan", "tyrosine"],
                               [0.6, 0.5, 0.5],
                               2)

        with self.assertRaises(AssertionError):
            GeneralParasResult(ParasModel.PARAS,
                               "ILIKEDATAEVENFAKEDATADIDISAYILIKED",
                               ["tryptophan", "tyrosine", "phenylalanine"],
                               [0.6, 0.5, 0.5],
                               5)

    def test_html(self) -> None:
        self.paras_result.uncertain = False
        html = self.paras_result.as_html()
        assert isinstance(html, Markup)
        assert "low confidence prediction" not in str(html)

        self.paras_result.uncertain = True
        html = self.paras_result.as_html()
        assert "low confidence prediction" in str(html)

    def test_str(self) -> None:
        assert "ParasResult" in str(self.paras_result)
        assert "ParasectResult" in str(self.parasect_result)

    def assertNearlyEqual(self, list_1, list_2):
        for i, element_1 in enumerate(list_1):
            element_2 = list_2[i]
            if not isclose(element_1, element_2, abs_tol=0.00001):
                self.fail(
                    f"Lists are not equal: {list_1}, {list_2}. \n First mismatching element: {i} ([{element_1}], [{element_2}])")


class TestParasResult(unittest.TestCase):
    def setUp(self):
        self.minimal_result = MockMinimalResult("ILIKEDATAEVENFAKEDATADIDISAYILIKED",
                                                ["tryptophan", "tyrosine", "phenylalanine"],
                                                [0.9, 0.05, 0.05])

        self.json_result = {"aa34": "ILIKEDATAEVENFAKEDATADIDISAYILIKED",
                            "top_substrates": ["tryptophan", "tyrosine"],
                            "top_confidences": [0.9, 0.05]}


    def test_from_paras_result(self):
        pred = ParasResult.from_paras_result(self.minimal_result, substrates_to_report=2)
        assert pred.aa34 == "ILIKEDATAEVENFAKEDATADIDISAYILIKED"
        assert pred.method == "paras"
        assert pred.predicted_substrate == "tryptophan"
        assert isclose(pred.confidence, 0.9)
        assert pred.top_substrates == ["tryptophan", "tyrosine"]
        assert len(pred.top_confidences) == 2

    def test_from_json(self):
        pred = ParasResult.from_json(self.json_result)
        assert pred.aa34 == "ILIKEDATAEVENFAKEDATADIDISAYILIKED"
        assert pred.method == "paras"
        assert pred.predicted_substrate == "tryptophan"
        assert isclose(pred.confidence, 0.9)
        assert pred.top_substrates == ["tryptophan", "tyrosine"]
        assert len(pred.top_confidences) == 2


class TestParasectResult(unittest.TestCase):
    def setUp(self):
        self.minimal_result = MockMinimalResult("ILIKEDATAEVENFAKEDATADIDISAYILIKED",
                                                ["tryptophan", "tyrosine", "phenylalanine"],
                                                [0.9, 0.8, 0.7])

        self.json_result = {"aa34": "ILIKEDATAEVENFAKEDATADIDISAYILIKED",
                            "top_substrates": ["tryptophan", "tyrosine"],
                            "top_confidences": [0.9, 0.8]}


    def test_from_parasect_result(self):
        pred = ParasectResult.from_parasect_result(self.minimal_result, substrates_to_report=2)
        assert pred.aa34 == "ILIKEDATAEVENFAKEDATADIDISAYILIKED"
        assert pred.method == "parasect"
        assert pred.predicted_substrate == "tryptophan"
        assert isclose(pred.confidence, 0.9)
        assert pred.top_substrates == ["tryptophan", "tyrosine"]
        assert len(pred.top_confidences) == 2

    def test_from_json(self):
        pred = ParasectResult.from_json(self.json_result)
        assert pred.aa34 == "ILIKEDATAEVENFAKEDATADIDISAYILIKED"
        assert pred.method == "parasect"
        assert pred.predicted_substrate == "tryptophan"
        assert isclose(pred.confidence, 0.9)
        assert pred.top_substrates == ["tryptophan", "tyrosine"]
        assert len(pred.top_confidences) == 2
