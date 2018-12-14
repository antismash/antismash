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
        a_domains = [("bpsA", 1), ("bpsA", 2), ("bpsA", 3),
                     ("bpsB", 1), ("bpsB", 2), ("bpsB", 3),
                     ("bpsC", 1),
                     ("bpsD", 1)]
        nrps_names = ['nrpspksdomains_%s_AMP-binding.%d' % a_dom for a_dom in a_domains]
        feature_names = nrps_names + ["nrpspksdomains_pks_CAL_domain.1"]
        assert set(results.domain_predictions) == set(feature_names)

        assert set(results.domain_predictions[feature_names[0]]) == {"NRPSPredictor2"}
        nrpspred2_results = {}
        for domain, methods in results.domain_predictions.items():
            if "CAL" in domain:
                continue
            nrpspred2_results[domain] = methods["NRPSPredictor2"].get_classification()
        expected_preds = [["leu"], ["bht"], ["asn"], ["hpg"], ["hpg"], ["bht"], ["dhpg"], ["tyr"], ["pk"]]
        expected_nrps2 = {name: pred for name, pred in zip(nrps_names, expected_preds)}
        assert nrpspred2_results == expected_nrps2

        cal = results.domain_predictions["nrpspksdomains_pks_CAL_domain.1"]["minowa_cal"]
        assert len(cal.predictions) == 5
        assert cal.predictions[0] == ["AHBA", 167.0]


        assert len(results.region_predictions[1]) == 2
        # as does this, though it still won't use domain docking
        pred = results.region_predictions[1][0]
        monomers = '(leu - bht - asn) + (hpg - hpg - bht) + (dhpg) + (tyr) + (pk)'
        assert pred.polymer == monomers
        assert not pred.domain_docking_used

        pred = results.region_predictions[1][1]
        assert pred.polymer == "(tyr) + (pk)"
        assert not pred.domain_docking_used

    def test_cp002271_c19(self):
        filename = path.get_full_path(__file__, 'data', 'CP002271.1.cluster019.gbk')
        results = helpers.run_and_regenerate_results_for_module(filename, nrps_pks, self.options)
        # catch ordering changes along with ensuring ATResults are there
        pred = results.domain_predictions["nrpspksdomains_STAUR_3982_PKS_AT.1"]
        assert pred["signature"].predictions[0][1].score == 87.5
        # ensure all genes are present and have the right consensus
        assert results.consensus == {'nrpspksdomains_STAUR_3982_PKS_AT.1': 'ohmmal',
                                     'nrpspksdomains_STAUR_3983_PKS_AT.1': 'ccmmal',
                                     'nrpspksdomains_STAUR_3984_PKS_AT.1': 'ccmmal',
                                     'nrpspksdomains_STAUR_3985_PKS_AT.1': 'pk',
                                     'nrpspksdomains_STAUR_3985_PKS_AT.2': 'ccmmal'}
        assert len(results.region_predictions) == 1
        assert list(results.region_predictions) == [1]
        assert len(results.region_predictions[1]) == 1
        # check the gene ordering and, in this case, that it used domain docking
        sc_pred = results.region_predictions[1][0]
        assert sc_pred.polymer == '(ccmmal) + (ccmmal) + (pk - ccmmal) + (ohmmal)'
        assert sc_pred.domain_docking_used
        assert len(results.domain_predictions) == 10
        expected_domains = {'nrpspksdomains_STAUR_3982_PKS_AT.1',
                            'nrpspksdomains_STAUR_3983_PKS_AT.1',
                            'nrpspksdomains_STAUR_3984_PKS_AT.1',
                            'nrpspksdomains_STAUR_3985_PKS_AT.1',
                            'nrpspksdomains_STAUR_3985_PKS_AT.2',
                            'nrpspksdomains_STAUR_3972_PKS_KR.1',
                            'nrpspksdomains_STAUR_3984_PKS_KR.1',
                            'nrpspksdomains_STAUR_3985_PKS_KR.1',
                            'nrpspksdomains_STAUR_3983_PKS_KR.1',
                            'nrpspksdomains_STAUR_3983_PKS_KR.1',
                            'nrpspksdomains_STAUR_3982_PKS_KR.1'}
        assert set(results.domain_predictions) == expected_domains
