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
        self.options = build_config(["--minimal"], isolated=True, modules=antismash.get_all_modules())
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
        results = helpers.run_and_regenerate_results_for_module(filename, nrps_pks, self.options)
        assert len(results.domain_predictions) == 9
        a_domains = ["bpsA_A1", "bpsA_A2", "bpsA_A3", "bpsB_A1", "bpsB_A2", "bpsB_A3",
                     "bpsC_A1", "bpsD_A1"]
        feature_names = ['nrpspksdomains_%s' % a_dom for a_dom in a_domains]
        feature_names.append("nrpspksdomains_pks_CAL1")
        assert set(results.domain_predictions) == set(feature_names)

        assert set(results.domain_predictions[feature_names[0]]) == {"NRPSPredictor2"}
        nrpspred2_results = {}
        for domain, methods in results.domain_predictions.items():
            if "CAL" in domain:
                continue
            nrpspred2_results[domain] = methods["NRPSPredictor2"].get_classification()
        expected_preds = [["leu"], ["bht"], ["asn"], ["hpg"], ["hpg"], ["bht"], ["dhpg"], ["tyr"], ["pk"]]
        expected_nrps2 = {name: pred for name, pred in zip(feature_names, expected_preds) if name[-2] == "A"}
        assert nrpspred2_results == expected_nrps2

        cal = results.domain_predictions["nrpspksdomains_pks_CAL1"]["minowa_cal"]
        assert len(cal.predictions) == 5
        assert cal.predictions[0] == ["AHBA", 167.0]

        # as does this, though it still won't use domain docking
        monomers = '(leu-bht-asn) + (hpg-hpg-bht) + (dhpg) + (tyr) + (pk)'
        assert results.cluster_predictions == {1: [monomers, False]}

    def test_cp002271_c19(self):
        filename = path.get_full_path(__file__, 'data', 'CP002271.1.cluster019.gbk')
        results = helpers.run_and_regenerate_results_for_module(filename, nrps_pks, self.options)
        # catch ordering changes along with ensuring ATResults are there
        assert results.domain_predictions["nrpspksdomains_STAUR_3982_AT1"]["signature"].predictions[0].score == 87.5
        # ensure all genes are present and have the right consensus
        assert results.consensus == {'nrpspksdomains_STAUR_3982_AT1': 'ohmmal',
                                     'nrpspksdomains_STAUR_3983_AT1': 'ccmmal',
                                     'nrpspksdomains_STAUR_3984_AT1': 'ccmmal',
                                     'nrpspksdomains_STAUR_3985_AT1': 'mmal',
                                     'nrpspksdomains_STAUR_3985_AT2': 'pk'}
        # check the gene ordering and, in this case, that it used domain docking
        assert results.cluster_predictions == {1: [
                '(ccmmal) + (ccmmal) + (mmal-pk) + (ohmmal)', True]}
        assert len(results.domain_predictions) == 10
        expected_domains = {'nrpspksdomains_STAUR_3982_AT1',
                            'nrpspksdomains_STAUR_3983_AT1',
                            'nrpspksdomains_STAUR_3984_AT1',
                            'nrpspksdomains_STAUR_3985_AT1',
                            'nrpspksdomains_STAUR_3985_AT2',
                            'nrpspksdomains_STAUR_3972_KR1',
                            'nrpspksdomains_STAUR_3984_KR1',
                            'nrpspksdomains_STAUR_3985_KR1',
                            'nrpspksdomains_STAUR_3983_KR1',
                            'nrpspksdomains_STAUR_3983_KR1',
                            'nrpspksdomains_STAUR_3982_KR1'}
        assert set(results.domain_predictions) == expected_domains
