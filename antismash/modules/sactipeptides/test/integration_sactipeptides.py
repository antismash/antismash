# License: GNU Affero General Public License v3 or later
# A copy of GNU AGPL v3 should have been included in this software package in LICENSE.txt.

# for test files, silence irrelevant and noisy pylint warnings
# pylint: disable=no-self-use,protected-access,missing-docstring

import unittest

import antismash
from antismash.common import path
from antismash.common.test import helpers
from antismash.config import build_config, update_config, destroy_config
from antismash.modules import sactipeptides
from antismash.modules.sactipeptides import SactiResults


class IntegrationSactipeptides(unittest.TestCase):
    def setUp(self):
        self.options = build_config(["--minimal", "--enable-sactipeptides"],
                                    isolated=True, modules=antismash.get_all_modules())

    def tearDown(self):
        destroy_config()

    def set_fimo_enabled(self, val):
        update_config({"without_fimo": not val})

    def gather_all_motifs(self, result):
        motifs = []
        for locii in result.clusters.values():
            for locus in locii:
                motifs.extend(result.motifs_by_locus[locus])
        return motifs

    def test_ap012495_end_to_end(self):
        result = helpers.run_and_regenerate_results_for_module(path.get_full_path(__file__, "data", "AP012495.1_c14.gbk"),
                        sactipeptides, self.options, expected_record_count=1)
        assert isinstance(result, SactiResults)
        assert list(result.motifs_by_locus) == ["BEST7613_6887"]
        prepeptide = result.motifs_by_locus["BEST7613_6887"][0]
        assert prepeptide.location.start == 9735
        assert prepeptide.location.end == 9867
        assert prepeptide.get_name() == "BEST7613_6887"
        assert prepeptide.leader == "MKKAVIVENK"
        assert prepeptide.core == "GCATCSIGAACLVDGPIPDFEIAGATGLFGLWG"
        self.assertAlmostEqual(prepeptide.score, 31.)

    def test_ap012495_end_to_end_all_orfs(self):
        # make sure that unannotated orfs are found if they are the precursor
        result = helpers.run_and_regenerate_results_for_module(path.get_full_path(__file__,
                                             "data", "AP012495.1_c14_missing_precursor.gbk"),
                        sactipeptides, self.options, expected_record_count=1)
        assert isinstance(result, SactiResults)
        assert list(result.motifs_by_locus) == ["allorf041"]
        prepeptide = result.motifs_by_locus["allorf041"][0]
        assert prepeptide.location.start == 9735
        assert prepeptide.location.end == 9867
        assert prepeptide.get_name() == "allorf041"
        assert prepeptide.leader == "MKKAVIVENK"
        assert prepeptide.core == "GCATCSIGAACLVDGPIPDFEIAGATGLFGLWG"
        self.assertAlmostEqual(prepeptide.score, 28.)
