# License: GNU Affero General Public License v3 or later
# A copy of GNU AGPL v3 should have been included in this software package in LICENSE.txt.

# for test files, silence irrelevant and noisy pylint warnings
# pylint: disable=no-self-use,protected-access,missing-docstring

import os
from tempfile import TemporaryDirectory
import unittest

import antismash
from antismash.main import detect_signature_genes, read_data, regenerate_results_for_record
from antismash.common.module_results import ModuleResults
from antismash.common.record_processing import parse_input_sequence
import antismash.common.test.helpers as helpers
from antismash.config import get_config, update_config, destroy_config, build_config
from antismash.modules import tta


class TtaIntegrationTest(unittest.TestCase):
    def setUp(self):
        options = build_config(["--minimal", "--tta"], isolated=True,
                                modules=antismash.get_all_modules())
        self.old_config = get_config().__dict__
        self.options = update_config(options)

    def tearDown(self):
        destroy_config()
        update_config(self.old_config)

    def test_nisin(self):
        record = parse_input_sequence(helpers.get_path_to_nisin_genbank())[0]
        detect_signature_genes(record, self.options)
        clusters = record.get_clusters()
        assert clusters
        for cluster in clusters:
            assert cluster.cds_children
        assert record.get_cds_features_within_clusters()
        assert record.get_feature_count() == 26

        assert tta.check_prereqs() == []
        assert tta.check_options(self.options) == []
        assert tta.is_enabled(self.options)
        prior_results = None
        results = tta.run_on_record(record, prior_results, self.options)
        assert isinstance(results, ModuleResults)
        assert len(results.features) == 174
        assert record.get_feature_count() == 26
        results.add_to_record(record)
        assert record.get_feature_count() == 200

    def test_nisin_complete(self):
        with TemporaryDirectory() as output_dir:
            args = ["--minimal", "--tta", "--output-dir", output_dir, helpers.get_path_to_nisin_genbank()]
            options = build_config(args, isolated=True, modules=antismash.get_all_modules())
            antismash.run_antismash(helpers.get_path_to_nisin_genbank(), options)

            # regen the results
            update_config({"reuse_results": os.path.join(output_dir, "nisin.json")})
            prior_results = read_data(None, options)
            record = prior_results.records[0]
            results = prior_results.results[0]
            regenned = regenerate_results_for_record(record, options, [tta], results)
            tta_results = regenned["antismash.modules.tta"]
            assert isinstance(tta_results, tta.TTAResults)
            assert len(tta_results.features) == 174
