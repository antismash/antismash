# License: GNU Affero General Public License v3 or later
# A copy of GNU AGPL v3 should have been included in this software package in LICENSE.txt.

import unittest

from antismash.main import get_all_modules, detect_signature_genes
from antismash.common import deprecated
from antismash.common.module_results import ModuleResults
import antismash.common.test.helpers as helpers
from antismash.config import args, get_config, update_config
from antismash.modules import tta

class TtaIntegrationTest(unittest.TestCase):
    def setUp(self):
        options = args.build_parser(modules=get_all_modules()).parse_args(["--tta"])
        self.old_config = get_config().__dict__
        self.options = update_config(options)

    def tearDown(self):
        update_config({})

    def test_nisin(self):
        record = deprecated.parse_input_sequence(helpers.get_path_to_nisin_genbank(), self.options)[0]
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
