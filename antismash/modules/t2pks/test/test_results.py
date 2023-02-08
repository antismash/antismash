# License: GNU Affero General Public License v3 or later
# A copy of GNU AGPL v3 should have been included in this software package in LICENSE.txt.

# for test files, silence irrelevant and noisy pylint warnings
# pylint: disable=use-implicit-booleaness-not-comparison,protected-access,missing-docstring

import json
import unittest

from antismash.modules.t2pks.results import (
    CDSPrediction,
    ProtoclusterPrediction,
    Prediction,
    T2PKSResults
)


class TestPrediction(unittest.TestCase):
    def test_stringification(self):
        pred = Prediction("somename", 15., 1e-89)
        assert str(pred) == "somename (Score: 15.0; E-value: 1e-89)"

    def test_reconstruction(self):
        pred = Prediction("somename", 15., 1e-89)
        reconstructed = Prediction.from_json(json.loads(json.dumps(pred.to_json())))
        assert pred is not reconstructed
        assert pred == reconstructed

    def test_invalid(self):
        with self.assertRaisesRegex(ValueError, "could not convert string to float"):
            Prediction("name", "not a float", 1e-89)
        with self.assertRaisesRegex(ValueError, "could not convert string to float"):
            Prediction("name", 1e-89, "not a float")

    def test_equality(self):
        first = Prediction("somename", 15., 1e-89)
        second = Prediction("somename", 15., 1e-89)
        assert first == second

        second = Prediction("tomename", 15., 1e-89)
        assert first != second

        second = Prediction("somename", 1., 1e-89)
        assert first != second

        second = Prediction("tomename", 15., 1e-88)
        assert first != second

        assert first != "not a Prediction"


class TestCDSPrediction(unittest.TestCase):
    def test_stringification_without_function(self):
        pred = CDSPrediction(protein_type="HAL", protein_function=None,
                             bitscore=20.5, evalue=1e-16)
        assert str(pred) == "HAL (Score: 20.5; E-value: 1e-16)"

    def test_stringification_with_function(self):
        pred = CDSPrediction(protein_type="CYC", protein_function="C7-C12",
                             bitscore=20.5, evalue=1e-16)
        assert str(pred) == "CYC C7-C12 (Score: 20.5; E-value: 1e-16)"

    def test_reconstruction(self):
        pred = CDSPrediction(protein_type="CYC", protein_function="C7-C12",
                             bitscore=20.5, evalue=1e-16)
        reconstructed = CDSPrediction.from_json(json.loads(json.dumps(pred.to_json())))
        assert str(reconstructed) == str(pred)

    def test_invalid_reconstruction(self):
        good = ["CYC", "C7-C12", 20.5, 1e-16]
        # ensure it works
        CDSPrediction.from_json(good)
        for i in range(4):
            bad = good[:i] + good[i+1:]
            # ensure it doesn't
            with self.assertRaisesRegex(ValueError, "Invalid CDSPrediction JSON"):
                CDSPrediction.from_json(bad)


def build_dummy_cds_predictions():
    preds = [CDSPrediction(protein_type="CYC", protein_function="C7-C12",
                           bitscore=20.5, evalue=1e-16),
             CDSPrediction(protein_type="HAL", protein_function=None,
                           bitscore=21.0, evalue=1e-160)]
    return {"cds_name": preds}


class TestProtoclusterPrediction(unittest.TestCase):
    def test_reconstruction(self):
        preds_by_cds = build_dummy_cds_predictions()

        starters = [Prediction('acetyl', 0., 0.)]
        elongations = [Prediction('7', 743.5, 1.2e-226)]
        classes = {"benzoisochromanequinone"}
        weights = {"acetyl_7": 451.23}
        cluster_pred = ProtoclusterPrediction(preds_by_cds,
                                              starter_units=starters,
                                              malonyl_elongations=elongations,
                                              product_classes=classes,
                                              molecular_weights=weights,
                                              start=100,
                                              end=2000
                                              )
        assert cluster_pred.cds_predictions == preds_by_cds
        assert cluster_pred.starter_units == starters
        assert cluster_pred.malonyl_elongations == elongations
        assert cluster_pred.product_classes == classes
        assert cluster_pred.molecular_weights == weights
        assert cluster_pred.start == 100
        assert cluster_pred.end == 2000
        reconstructed = ProtoclusterPrediction.from_json(json.loads(json.dumps(cluster_pred.to_json())))
        assert cluster_pred.cds_predictions == preds_by_cds
        assert reconstructed.starter_units == starters
        assert reconstructed.malonyl_elongations == elongations
        assert reconstructed.product_classes == classes
        assert reconstructed.molecular_weights == weights
        assert reconstructed.start == 100
        assert reconstructed.end == 2000


def build_dummy_cluster_prediction():
    return ProtoclusterPrediction(build_dummy_cds_predictions(),
                                  [Prediction('acetyl', 0., 0.)],
                                  [Prediction('7', 743.5, 1.2e-226)],
                                  {"benzoisochromanequinone"},
                                  {"acetyl_7": 451.23},
                                  start=0, end=1500)


class TestT2PKSResults(unittest.TestCase):
    def test_reconstruction(self):
        results = T2PKSResults("rec_id")
        results.cluster_predictions[7] = build_dummy_cluster_prediction()
        reconstructed = T2PKSResults.from_json(json.loads(json.dumps(results.to_json())), _record=None)
        assert str(results) == str(reconstructed)

        assert list(reconstructed.cluster_predictions) == [7]
        assert str(reconstructed.cluster_predictions[7]) == str(results.cluster_predictions[7])

        original_score = results.cluster_predictions[7].cds_predictions["cds_name"][0].bitscore
        recon_score = reconstructed.cluster_predictions[7].cds_predictions["cds_name"][0].bitscore
        assert original_score == recon_score

    def test_mismatching_schema(self):
        results = T2PKSResults("rec_id")
        results.cluster_predictions[7] = build_dummy_cluster_prediction()
        results_json = json.loads(json.dumps(results.to_json()))
        results_json["schema_version"] += 1
        reconstructed = T2PKSResults.from_json(results_json, _record=None)
        assert reconstructed is None
