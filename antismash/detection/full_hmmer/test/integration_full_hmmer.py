# License: GNU Affero General Public License v3 or later
# A copy of GNU AGPL v3 should have been included in this software package in LICENSE.txt.

# for test files, silence irrelevant and noisy pylint warnings
# pylint: disable=no-self-use,protected-access,missing-docstring

import unittest

import antismash
from antismash.common import record_processing
from antismash.common.test import helpers
from antismash.config import build_config, destroy_config, update_config
from antismash.detection import full_hmmer


class TestFullHmmer(unittest.TestCase):
    def setUp(self):
        self.options = build_config(["--fh-analysis", "--minimal",
                                     "--fh-min-score", "1",
                                     "--fh-max-evalue", "0.01"],
                                    isolated=True,
                                    modules=antismash.get_all_modules())

    def tearDown(self):
        destroy_config()

    def check_add_to_record(self, input_file, results):
        record = record_processing.parse_input_sequence(input_file)[0]
        assert not record.get_pfam_domains()
        results.add_to_record(record)
        assert len(record.get_pfam_domains()) == len(results.hits)

    def test_reuse(self):
        nisin = helpers.get_path_to_nisin_genbank()
        record = record_processing.parse_input_sequence(nisin)[0]

        results = helpers.run_and_regenerate_results_for_module(nisin, full_hmmer,
                                         self.options, expected_record_count=1)
        json = results.to_json()
        assert len(results.hits) == 24
        self.check_add_to_record(nisin, results)

        # test regeneration when thresholds are less restrictive
        original_score_threshold = self.options.fh_min_score
        original_evalue_threshold = self.options.fh_max_evalue

        new_score_threshold = original_score_threshold / 2
        update_config({"fh_min_score": new_score_threshold})
        new_results = full_hmmer.full_hmmer.FullHmmerResults.from_json(json, record)
        assert new_results is None
        update_config({"fh_min_score": original_score_threshold})

        new_evalue_threshold = original_evalue_threshold * 2
        update_config({"fh_max_evalue": new_evalue_threshold})
        new_results = full_hmmer.full_hmmer.FullHmmerResults.from_json(json, record)
        assert new_results is None
        update_config({"fh_max_evalue": original_evalue_threshold})

        # test regeneration when evalue threshold is more restrictive
        new_evalue_threshold = sorted(hit["evalue"] for hit in results.hits)[12]
        assert new_evalue_threshold < self.options.fh_max_evalue
        new_hits = []
        for hit in results.hits:
            if hit["evalue"] <= new_evalue_threshold:
                new_hits.append(hit)
        new_hits.sort(key=lambda x: x["evalue"])
        assert len(new_hits) < 24

        update_config({"fh_max_evalue": new_evalue_threshold})
        new_results = full_hmmer.full_hmmer.FullHmmerResults.from_json(json, record)
        update_config({"fh_max_evalue": original_evalue_threshold})
        assert sorted(new_results.hits, key=lambda x: x["evalue"]) == new_hits
        self.check_add_to_record(nisin, results)

        # test regeneration when score threshold is more restrictive
        new_score_threshold = sorted(hit["score"] for hit in results.hits)[12]
        assert new_score_threshold > self.options.fh_min_score
        new_hits = []
        for hit in results.hits:
            if hit["score"] >= new_score_threshold:
                new_hits.append(hit)
        new_hits.sort(key=lambda x: x["score"])
        assert len(new_hits) < 24

        update_config({"fh_min_score": new_score_threshold})
        new_results = full_hmmer.full_hmmer.FullHmmerResults.from_json(json, record)
        update_config({"fh_min_score": original_score_threshold})
        assert sorted(new_results.hits, key=lambda x: x["score"]) == new_hits
        self.check_add_to_record(nisin, results)
