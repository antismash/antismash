# License: GNU Affero General Public License v3 or later
# A copy of GNU AGPL v3 should have been included in this software package in LICENSE.txt.

import unittest

from antismash.main import gather_modules, detect_signature_genes
from antismash.common import deprecated
from antismash.common.module_results import ModuleResults
import antismash.common.test.helpers as helpers
from antismash.config import args
from antismash.modules import tta

class TtaIntegrationTest(unittest.TestCase):
    def setUp(self):
        options = args.build_parser(modules=gather_modules(with_genefinding=True)).parse_args(["--tta"])
        self.old_config = args.Config().__dict__
        self.options = args.Config(options)

    def tearDown(self):
        args.Config({})

    def test_nisin(self):
        record = deprecated.parse_input_sequence(helpers.get_path_to_nisin_genbank(), self.options)[0]
        detect_signature_genes(record, self.options)
        clusters = record.get_clusters()
        assert clusters
        for cluster in clusters:
            assert cluster.cds_children
        assert deprecated.get_cds_features_within_clusters(record)

        assert tta.check_prereqs() == []
        assert tta.check_options(self.options) == []
        assert tta.is_enabled(self.options)
        results = tta.run_on_record(record, self.options)
        assert isinstance(results, ModuleResults)
        assert len(results.features) == 174
