# License: GNU Affero General Public License v3 or later
# A copy of GNU AGPL v3 should have been included in this software package in LICENSE.txt.

# for test files, silence irrelevant and noisy pylint warnings
# pylint: disable=use-implicit-booleaness-not-comparison,protected-access,missing-docstring

import json
import unittest
from unittest.mock import patch

from antismash.main import get_all_modules
from antismash.common import path, subprocessing
from antismash.common.secmet import Record
from antismash.common.test import helpers
from antismash.config import build_config, destroy_config, update_config
from antismash.modules import cluster_compare


DUMMY_MIBIG_CONFIG = path.get_full_path(__file__, "data", "dummy_mibig.json")
DUMMY_MIBIG_LOCATION = path.get_full_path(__file__, "data", "dummy_mibig_db")


class TestFull(unittest.TestCase):
    def setUp(self):
        self.options = build_config([], isolated=True, modules=get_all_modules())

    def tearDown(self):
        destroy_config()

    def test_result_conversion(self):
        update_config({"cc_mibig": True})
        nisin = Record.from_genbank(helpers.get_path_to_nisin_with_detection())[0]
        with open(path.get_full_path(__file__, "data", "nisin.out"), encoding="utf-8") as handle:
            trimmed_output = handle.read()
        with patch.object(subprocessing, "run_diamond_search", return_value=trimmed_output):
            results = cluster_compare.run_on_record(nisin, None, self.options)
        assert results.by_database["MIBiG"].by_region[1]
        # ensure JSON conversion of results gives the same result
        raw = json.loads(json.dumps(results.to_json()))
        regenerated = cluster_compare.regenerate_previous_results(raw, nisin, self.options)
        regen_raw = json.loads(json.dumps(regenerated.to_json()))
        assert regen_raw == raw

    def test_full(self):
        update_config({"cc_custom_dbs": [DUMMY_MIBIG_CONFIG]})
        assert cluster_compare.is_enabled(self.options)
        input_path = helpers.get_path_to_balhymicin_genbank()
        config = cluster_compare.DBConfig.from_file(DUMMY_MIBIG_CONFIG, DUMMY_MIBIG_LOCATION)
        with patch.object(cluster_compare.DBConfig, "from_json", return_value=config):
            results = helpers.run_and_regenerate_results_for_module(input_path, cluster_compare, self.options)

        proto = results.by_database["MIBiG"].by_region[1]["ProtoToRegion_RiQ"]
        # ordering should be correct
        assert proto.scores_by_region[0][1] > proto.scores_by_region[1][1]
        # check the winner makes sense for proto to region, the value won't be 1.0
        assert proto.scores_by_region[0][0].accession == "BGC0000311"
        self.assertAlmostEqual(proto.scores_by_region[0][1], 1.80874, places=5)

        region = results.by_database["MIBiG"].by_region[1]["RegionToRegion_RiQ"]
        # ordering should be correct
        assert region.scores_by_region[0][1] > region.scores_by_region[1][1]
        # again, check the winner, but this time it *should* be 1.0 as there's no protoclusters used
        assert region.scores_by_region[0][0].accession == "BGC0000311", region.scores_by_region
        self.assertAlmostEqual(region.scores_by_region[0][1], 1.0, places=5)
