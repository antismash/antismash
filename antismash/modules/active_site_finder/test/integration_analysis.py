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
from antismash.modules import active_site_finder


class TestAnalyses(unittest.TestCase):
    def setUp(self):
        # skipping clusterhmmer and the p450 potential hits for speed
        self.options = build_config(["--asf", "--minimal", "--enable-html"],
                                    isolated=True,
                                    modules=antismash.get_all_modules())

    def test_regeneration(self):
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
