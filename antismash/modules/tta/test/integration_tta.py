# License: GNU Affero General Public License v3 or later
# A copy of GNU AGPL v3 should have been included in this software package in LICENSE.txt.

# for test files, silence irrelevant and noisy pylint warnings
# pylint: disable=use-implicit-booleaness-not-comparison,protected-access,missing-docstring

import os
from tempfile import TemporaryDirectory
import unittest

import antismash
from antismash.main import read_data
from antismash.common.module_results import ModuleResults
from antismash.common.secmet import Record
from antismash.common.test import helpers
from antismash.config import get_config, update_config, destroy_config, build_config
from antismash.modules import tta


class TtaIntegrationTest(unittest.TestCase):
    def setUp(self):
        options = build_config(["--minimal", "--enable-tta", "--tta-threshold", "0"],
                               isolated=True, modules=antismash.get_all_modules())
        self.old_config = get_config().__dict__
        self.options = update_config(options)

    def tearDown(self):
        destroy_config()
        update_config(self.old_config)

    def test_nisin(self):
        record = Record.from_genbank(helpers.get_path_to_nisin_with_detection(), taxon="bacteria")[0]
        regions = record.get_regions()
        assert regions
        for region in regions:
            assert region.cds_children
        assert record.get_cds_features_within_regions()
        before_count = record.get_feature_count()

        assert tta.check_prereqs(self.options) == []
        assert tta.check_options(self.options) == []
        assert tta.is_enabled(self.options)
        prior_results = None
        results = tta.run_on_record(record, prior_results, self.options)
        assert isinstance(results, ModuleResults)
        assert len(results.features) == 174
        assert record.get_feature_count() == before_count
        results.add_to_record(record)
        assert record.get_feature_count() == before_count + 174

    def test_nisin_complete(self):
        with TemporaryDirectory() as output_dir:
            args = ["--minimal", "--enable-tta", "--tta-threshold", "0",
                    "--output-dir", output_dir, helpers.get_path_to_nisin_genbank()]
            options = build_config(args, isolated=True, modules=antismash.get_all_modules())
            antismash.run_antismash(helpers.get_path_to_nisin_genbank(), options)

            # regen the results
            update_config({"reuse_results": os.path.join(output_dir, "nisin.json")})
            prior_results = read_data(None, options)
            record = prior_results.records[0]
            results = prior_results.results[0]
            tta_results = tta.regenerate_previous_results(results.get("antismash.modules.tta"), record, options)
            assert isinstance(tta_results, tta.TTAResults)
            assert len(tta_results.features) == 174

            # raise the threshold above the gc_content and ensure regenned has no hits
            update_config({"tta_threshold": 0.65})
            tta_results = tta.regenerate_previous_results(results.get("antismash.modules.tta"), record, options)
            assert isinstance(tta_results, tta.TTAResults)
            assert not tta_results.features
