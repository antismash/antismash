# License: GNU Affero General Public License v3 or later
# A copy of GNU AGPL v3 should have been included in this software package in LICENSE.txt.

# for test files, silence irrelevant and noisy pylint warnings
# pylint: disable=use-implicit-booleaness-not-comparison,protected-access,missing-docstring

import json
import unittest
from unittest.mock import patch

from antismash.common.secmet.test.helpers import DummyAntismashDomain
from antismash.modules.nrps_pks import paras


class TestParasResult(unittest.TestCase):
    def setUp(self) -> None:
        self.predictions = [
            (0.9, "alanine"),
            (0.8, "glycine"),
        ]

    def test_json_roundtrip(self):
        result = paras.ParasResult(self.predictions)
        data = json.dumps(result.to_json())
        restored = paras.ParasResult.from_json(json.loads(data))
        assert restored.get_classification() == result.get_classification()

    def test_html(self):
        no_preds = paras.ParasResult([])
        assert "No hits above threshold" in str(no_preds.as_html())

        preds = paras.ParasResult(self.predictions)
        html_output = str(preds.as_html())

        for score, name in self.predictions:
            assert name in html_output
            assert f"{score:.2f}" in html_output

    def test_classification(self):
        no_preds = paras.ParasResult([])
        assert [] == no_preds.get_classification()

        simple = paras.ParasResult(self.predictions)
        assert ["alanine"] == simple.get_classification()

        tie = paras.ParasResult([
            (0.9, "alanine"),
            (0.9, "glycine"),
            (0.8, "serine"),
        ])
        assert ["alanine", "glycine"] == tie.get_classification()


class TestParasRun(unittest.TestCase):
    def setUp(self) -> None:
        domain_a = DummyAntismashDomain(locus_tag="A", domain_id="A_1")
        domain_a._translation = "MAGICHAT"
        domain_b = DummyAntismashDomain(locus_tag="A", domain_id="A_2")
        domain_b._translation = "MAGICCAT"
        self.domains = [domain_a, domain_b]
        self.predictions = [
            [(0.9, "alanine"), (0.8, "glycine")],
            [(0.4, "serine"), (0.4, "threonine")],
        ]

    def test_run_paras(self):
        with patch.object(paras, "run_paras_bulk", return_value=self.predictions):
            predictions = paras.run_paras(self.domains)

        expected_results = {
            "A_1": ["alanine"],
            "A_2": ["serine", "threonine"],
        }
        for name, expected in expected_results.items():
            assert expected == predictions[name].get_classification()
