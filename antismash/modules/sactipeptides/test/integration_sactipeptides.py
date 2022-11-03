# License: GNU Affero General Public License v3 or later
# A copy of GNU AGPL v3 should have been included in this software package in LICENSE.txt.

# for test files, silence irrelevant and noisy pylint warnings
# pylint: disable=use-implicit-booleaness-not-comparison,protected-access,missing-docstring

import unittest

import antismash
from antismash.common import path
from antismash.common.test import helpers
from antismash.config import build_config, destroy_config
from antismash.modules import sactipeptides
from antismash.modules.sactipeptides import SactiResults


class IntegrationSactipeptides(unittest.TestCase):
    def setUp(self):
        self.options = build_config(["--minimal", "--enable-sactipeptides", "--enable-html"],
                                    isolated=True, modules=antismash.get_all_modules())

    def tearDown(self):
        destroy_config()

    def run_analyis(self, filename):
        data_file = path.get_full_path(__file__, "data", filename)
        return helpers.run_and_regenerate_results_for_module(data_file, sactipeptides, self.options)

    def test_ap012495_end_to_end(self):
        result = self.run_analyis("AP012495.1_c14.gbk")
        assert isinstance(result, SactiResults)
        assert list(result.motifs_by_locus) == ["BEST7613_6887"]
        prepeptide = result.motifs_by_locus["BEST7613_6887"][0]
        assert prepeptide.location.start == 9735
        assert prepeptide.location.end == 9867
        assert prepeptide.get_name() == "BEST7613_6887_sactipeptide"
        assert prepeptide.leader == "MKKAVIVENK"
        assert prepeptide.core == "GCATCSIGAACLVDGPIPDFEIAGATGLFGLWG"
        self.assertAlmostEqual(prepeptide.score, 33.)

    def test_ap012495_end_to_end_all_orfs(self):
        # make sure that unannotated orfs are found if they are the precursor
        result = self.run_analyis("AP012495.1_c14_missing_precursor.gbk")
        assert isinstance(result, SactiResults)
        assert list(result.motifs_by_locus) == ["allorf_09736_09867"]
        prepeptide = result.motifs_by_locus["allorf_09736_09867"][0]
        assert prepeptide.location.start == 9735
        assert prepeptide.location.end == 9867
        assert prepeptide.get_name() == "allorf_09736_09867_sactipeptide"
        assert prepeptide.leader == "MKKAVIVENK"
        assert prepeptide.core == "GCATCSIGAACLVDGPIPDFEIAGATGLFGLWG"
        self.assertAlmostEqual(prepeptide.score, 33.)
