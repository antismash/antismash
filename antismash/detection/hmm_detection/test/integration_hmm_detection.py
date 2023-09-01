# License: GNU Affero General Public License v3 or later
# A copy of GNU AGPL v3 should have been included in this software package in LICENSE.txt.

# for test files, silence irrelevant and noisy pylint warnings
# pylint: disable=use-implicit-booleaness-not-comparison,protected-access,missing-docstring

import json
import unittest

from antismash.common.secmet import Record
from antismash.common.secmet.test.helpers import DummySubRegion
from antismash.common.test.helpers import get_path_to_nisin_genbank
from antismash.config import build_config, destroy_config
from antismash.detection import hmm_detection


class TestSubregionAnnotations(unittest.TestCase):
    def setUp(self):
        self.options = build_config([], isolated=True, modules=[hmm_detection])

    def tearDown(self):
        destroy_config()
        hmm_detection._RULESETS.clear()

    def test_subregions_annotated(self):
        record = Record.from_genbank(get_path_to_nisin_genbank())[0]
        record.strip_antismash_annotations()
        assert not record.get_regions()
        assert not record.get_subregions()
         # this next hack requires that the ruleset cache is cleared in the test teardown
         # because it's nuking the default ruleset that all following tests (not just this file)
         # will use
        hmm_detection.get_ruleset(self.options)._rules = []

        results = hmm_detection.run_on_record(record, None, self.options)
        assert not results.get_predicted_protoclusters()
        for cds in ["nisB", "nisC"]:
            assert not record.get_cds_by_name(cds).sec_met
        cutoff = record.get_cds_by_name("nisB").location.end + 10
        record.add_subregion(DummySubRegion(end=cutoff))

        results = hmm_detection.run_on_record(record, None, self.options)
        assert record.get_cds_by_name("nisB").sec_met
        assert not record.get_cds_by_name("nisC").sec_met

        # and then a json conversion, even without subregion added
        record = Record.from_genbank(get_path_to_nisin_genbank())[0]
        record.strip_antismash_annotations()
        raw = json.loads(json.dumps(results.to_json()))
        hmm_detection.regenerate_previous_results(raw, record, self.options)
        assert record.get_cds_by_name("nisB").sec_met
        assert not record.get_cds_by_name("nisC").sec_met
