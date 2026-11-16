# License: GNU Affero General Public License v3 or later
# A copy of GNU AGPL v3 should have been included in this software package in LICENSE.txt.

# for test files, silence irrelevant and noisy pylint warnings
# pylint: disable=use-implicit-booleaness-not-comparison,protected-access,missing-docstring

import unittest

import antismash
from antismash.config import build_config
from antismash.common import secmet
from antismash.common.record_processing import parse_input_sequence
from antismash.common.test import helpers
from antismash.detection import nrps_pks_domains
from antismash.modules import active_site_finder
from antismash.outputs import html


class TestAnalyses(unittest.TestCase):
    def setUp(self):
        module_arg = helpers.build_module_limitation_arg([
            active_site_finder, html, nrps_pks_domains,
        ])
        self.options = build_config(module_arg + ["--asf"],
                                    isolated=True,
                                    modules=antismash.get_all_modules())

    def test_regeneration(self):
        assert active_site_finder.is_enabled(self.options)
        datafile = helpers.get_path_to_balhymicin_genbank()
        results = helpers.run_and_regenerate_results_for_module(datafile, active_site_finder, self.options)
        assert results.pairings
        for domain, labels in results.pairings:
            for label in labels:
                assert label
                assert isinstance(label, str)
            assert isinstance(domain, secmet.AntismashDomain)
        record = parse_input_sequence(datafile)

        # check the reuse portion works
        rerun = active_site_finder.run_on_record(record, results, self.options)
        assert rerun is results  # specifically checking it's the same object

        with self.assertRaisesRegex(AssertionError, "str"):
            active_site_finder.run_on_record(record, "invalid", self.options)
