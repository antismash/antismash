# License: GNU Affero General Public License v3 or later
# A copy of GNU AGPL v3 should have been included in this software package in LICENSE.txt.

# for test files, silence irrelevant and noisy pylint warnings
# pylint: disable=use-implicit-booleaness-not-comparison,protected-access,missing-docstring

import unittest
from unittest.mock import patch

from antismash.common.hmm_rule_parser.cluster_prediction import RuleDetectionResults
from antismash.common.hmm_rule_parser.structures import Multipliers
from antismash.common.secmet.test.helpers import DummyRecord
from antismash.detection.subclusters import results as results_module
from antismash.detection.subclusters.results import SubclusterDetectionResults


class TestFromJSON(unittest.TestCase):
    def setUp(self):
        self.record = DummyRecord(seq="A" * 100)
        self.json = {
            "schema_version": SubclusterDetectionResults.schema_version,
            "record_id": self.record.id,
            "strictness": "relaxed",
            "subregion_mode": "clip",
            "rule_names": ["b", "a"],
            "rule_results": {"dummy": True},
        }

    def test_schema_version(self):
        self.json["rule_results"] = RuleDetectionResults({}, "subclusters", [],
                                                         Multipliers()).to_json()
        assert SubclusterDetectionResults.from_json(self.json, self.record) is not None

        self.json["schema_version"] = -1
        assert SubclusterDetectionResults.from_json(self.json, self.record) is None

    def test_rule_results_discarded(self):
        with patch.object(results_module.RuleDetectionResults, "from_json",
                          return_value=None) as patched:
            assert SubclusterDetectionResults.from_json(self.json, self.record) is None
        patched.assert_called_once_with(self.json["rule_results"], self.record)
