# License: GNU Affero General Public License v3 or later
# A copy of GNU AGPL v3 should have been included in this software package in LICENSE.txt.

# for test files, silence irrelevant and noisy pylint warnings
# pylint: disable=use-implicit-booleaness-not-comparison,protected-access,missing-docstring

import unittest

import antismash
from antismash.common import record_processing, pfamdb
from antismash.common.test import helpers
from antismash.config import build_config, destroy_config
from antismash.detection import cluster_hmmer


class TestClusterHmmer(unittest.TestCase):
    def setUp(self):
        self.original_min_score = cluster_hmmer.MIN_SCORE
        self.original_max_evalue = cluster_hmmer.MAX_EVALUE
        self.options = build_config(["--clusterhmmer", "--minimal"],
                                    isolated=True,
                                    modules=antismash.get_all_modules())
        self.pfam_db = pfamdb.get_latest_db_path(self.options.database_dir)

    def tearDown(self):
        self.set_max_evalue(self.original_max_evalue)
        self.set_min_score(self.original_min_score)
        destroy_config()

    def set_max_evalue(self, evalue):
        cluster_hmmer.MAX_EVALUE = evalue

    def set_min_score(self, score):
        cluster_hmmer.MIN_SCORE = score

    def check_add_to_record(self, input_file, results):
        record = record_processing.parse_input_sequence(input_file)[0]
        assert not record.get_pfam_domains()
        results.add_to_record(record)
        assert len(record.get_pfam_domains()) == len(results.hits)

    def test_reuse(self):
        nisin = helpers.get_path_to_nisin_genbank()
        record = record_processing.parse_input_sequence(nisin)[0]

        results = helpers.run_and_regenerate_results_for_module(nisin, cluster_hmmer, self.options)
        json = results.to_json()
        assert len(results.hits) == 12
        self.check_add_to_record(nisin, results)

        # test regeneration when thresholds are less restrictive
        new_score_threshold = self.original_min_score - .1
        self.set_min_score(new_score_threshold)
        new_results = cluster_hmmer.regenerate_previous_results(json, record, self.options)
        assert new_results is None
        self.set_min_score(self.original_min_score)

        new_evalue_threshold = self.original_max_evalue + .1
        self.set_max_evalue(new_evalue_threshold)
        new_results = cluster_hmmer.regenerate_previous_results(json, record, self.options)
        assert new_results is None
        self.set_max_evalue(self.original_max_evalue)

        # test regeneration when evalue threshold is more restrictive
        new_evalue_threshold = sorted(hit.evalue for hit in results.hits)[6]
        assert new_evalue_threshold < self.original_max_evalue
        new_hits = []
        for hit in results.hits:
            if hit.evalue <= new_evalue_threshold:
                new_hits.append(hit)
        new_hits.sort(key=lambda x: x.evalue)
        assert len(new_hits) < 13

        self.set_max_evalue(new_evalue_threshold)
        new_results = cluster_hmmer.regenerate_previous_results(json, record, self.options)
        self.set_max_evalue(self.original_max_evalue)
        assert sorted(new_results.hits, key=lambda x: x.evalue) == new_hits
        self.check_add_to_record(nisin, results)

        # test regeneration when score threshold is more restrictive
        new_score_threshold = sorted(hit.score for hit in results.hits)[6]
        assert new_score_threshold > cluster_hmmer.MIN_SCORE
        new_hits = []
        for hit in results.hits:
            if hit.score >= new_score_threshold:
                new_hits.append(hit)
        new_hits.sort(key=lambda x: x.score)
        assert len(new_hits) < 13

        self.set_min_score(new_score_threshold)
        new_results = cluster_hmmer.regenerate_previous_results(json, record, self.options)
        self.set_min_score(self.original_min_score)
        assert sorted(new_results.hits, key=lambda x: x.score) == new_hits
        self.check_add_to_record(nisin, results)
