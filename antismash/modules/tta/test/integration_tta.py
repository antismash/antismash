# License: GNU Affero General Public License v3 or later
# A copy of GNU AGPL v3 should have been included in this software package in LICENSE.txt.

# for test files, silence irrelevant and noisy pylint warnings
# pylint: disable=no-self-use,protected-access,missing-docstring

import unittest

from antismash.main import get_all_modules, detect_signature_genes
from antismash.common.module_results import ModuleResults
from antismash.common.record_processing import parse_input_sequence
import antismash.common.test.helpers as helpers
from antismash.config import args, get_config, update_config, destroy_config
from antismash.modules import tta


class TtaIntegrationTest(unittest.TestCase):
    def setUp(self):
        options = args.build_parser(modules=get_all_modules()).parse_args(["--tta"])
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

        assert tta.check_prereqs() == []
        assert tta.check_options(self.options) == []
        assert tta.is_enabled(self.options)
        prior_results = None
        results = tta.run_on_record(record, prior_results, self.options)
        assert isinstance(results, ModuleResults)
        assert len(results.features) == 174
