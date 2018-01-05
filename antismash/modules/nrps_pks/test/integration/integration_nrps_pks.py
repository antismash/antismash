# License: GNU Affero General Public License v3 or later
# A copy of GNU AGPL v3 should have been included in this software package in LICENSE.txt.

# for test files, silence irrelevant and noisy pylint warnings
# pylint: disable=no-self-use,protected-access,missing-docstring

import unittest

from helperlibs.wrappers.io import TemporaryDirectory
import minimock

import antismash
from antismash.common import path
from antismash.common.test import helpers
from antismash.config import build_config, destroy_config
from antismash.modules import nrps_pks

class IntegrationWithoutNRPSPKS(unittest.TestCase):
    def setUp(self):
        self.options = build_config(["--minimal"],
                                    isolated=True, modules=antismash.get_all_modules())
        assert not nrps_pks.is_enabled(self.options)
        self.tracker = minimock.TraceTracker()
        minimock.mock("nrps_pks.run_on_record", tracker=self.tracker)

    def tearDown(self):
        destroy_config()
        minimock.restore()

    def test_minimal(self):
        with TemporaryDirectory(change=True) as tempdir:
            self.options = build_config(["--minimal", "--output-dir", tempdir],
                                    isolated=True, modules=antismash.get_all_modules())
            antismash.main.run_antismash(helpers.get_path_to_balhymicin_genbank(),
                                         self.options)
        # make sure it didn't run
        minimock.assert_same_trace(self.tracker, "")

class IntegrationNRPSPKS(unittest.TestCase):
    def setUp(self):
        self.options = build_config(["--minimal", "--enable-nrps-pks"],
                                    isolated=True, modules=antismash.get_all_modules())

    def tearDown(self):
        destroy_config()

    def test_balhymicin(self):
        filename = helpers.get_path_to_balhymicin_genbank()
        results = helpers.run_and_regenerate_results_for_module(filename, nrps_pks,
                                             self.options, expected_record_count=1)
        for key, val in results.pks.method_results.items():
            if key != "minowa_cal":
                assert not val
                continue
            assert len(val) == 1
            assert len(val["pks_CAL1"]) == 5
            assert val["pks_CAL1"][0] == ["AHBA", 167.0]
        # when the NRPS subsections are added, this needs to change
        assert results.nrps == {}
        # as does this, though it still won't use domain docking
        assert results.cluster_predictions == {'1': [
                '(nrp-nrp-nrp) + (nrp-nrp-nrp) + (nrp) + (nrp) + (pk)', False]}

    def test_CP002271_c19(self):
        filename = path.get_full_path(__file__, 'data', 'CP002271.1.cluster019.gbk')
        results = helpers.run_and_regenerate_results_for_module(filename, nrps_pks,
                                             self.options, expected_record_count=1)
        # catch ordering changes along with ensuring ATResults are there
        assert results.pks.method_results["signature"]["STAUR_3982_AT1"][0].score == 87.5
        # ensure all genes are present and have the right consensus
        assert results.consensus == {'STAUR_3982_AT1': 'ohmmal',
                                     'STAUR_3983_AT1': 'ccmmal',
                                     'STAUR_3984_AT1': 'ccmmal',
                                     'STAUR_3985_AT1': 'pk',
                                     'STAUR_3985_AT2': 'pk'}
        # check the gene ordering and, in this case, that it used domain docking
        assert results.cluster_predictions == {'1': [
                '(ccmmal) + (ccmmal) + (pk-pk) + (ohmmal)', True]}
        # no A domains in the cluster, so make sure no NRPS results
        assert results.nrps == {}
