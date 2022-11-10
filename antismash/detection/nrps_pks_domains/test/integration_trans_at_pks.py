# License: GNU Affero General Public License v3 or later
# A copy of GNU AGPL v3 should have been included in this software package in LICENSE.txt.

# for test files, silence irrelevant and noisy pylint warnings
# pylint: disable=no-self-use,protected-access,missing-docstring

import unittest

import antismash
from antismash.common.test import helpers
from antismash.config import build_config, destroy_config
from antismash.modules import nrps_pks
from antismash.detection import nrps_pks_domains


class IntegrationTRANSATPKS(unittest.TestCase):
    def setUp(self):
        self.options = build_config(["--minimal", "--enable-nrps-pks", "--enable-html"],
                                    isolated=True, modules=antismash.get_all_modules())

    def tearDown(self):
        destroy_config()

    def test_elansolid(self):
        filename = helpers.get_path_to_elansolid_genbank()
        results = helpers.run_and_regenerate_results_for_module(filename, nrps_pks, self.options)
        assert len(results.domain_predictions) == 14
        domains = [('elaJ_PKS_KR', 1), ('elaO_PKS_KR', 1),
                   ('elaJ_CAL_domain', 1), ('elaQ_PKS_KR', 1),
                   ('elaP_PKS_KR', 3), ('elaP_PKS_KR', 2),
                   ('elaK_PKS_KR', 2), ('elaK_PKS_KR', 3),
                   ('elaR_PKS_KR', 1), ('elaQ_PKS_KR', 2),
                   ('elaP_PKS_KR', 1), ('elaK_PKS_KR', 1),
                   ('elaB_PKS_AT', 1), ('elaC_PKS_AT', 1)]
        domain_names = ['nrpspksdomains_%s.%d' % dom for dom in domains]
        assert set(results.domain_predictions) == set(domain_names)
        assert set(results.domain_predictions[domain_names[0]]) == {"kr_activity"}
        pred = results.region_predictions[1][0]
        monomers = '(ccmal - Me-ohmal) + (ccmal - ccmal - mal) + (ohmal) + (ohmal - Me-ccmal) + (ohmal - Me-ccmal - ccmal) + (ohmal) + (ohmal - Me-ccmal) + (ohmal - Me-ccmal - ccmal)'
        assert pred.polymer == monomers
        results_domains = helpers.run_and_regenerate_results_for_module(filename, nrps_pks_domains, self.options)
        ks_subtypes = []
        ks_subtypes_true = ['Enediyne-KS', 'Trans-AT-KS', 'Trans-AT-KS', 'Trans-AT-KS',
                            'Trans-AT-KS', 'Trans-AT-KS', 'Trans-AT-KS', 'Trans-AT-KS',
                            'Trans-AT-KS', 'Trans-AT-KS', 'Trans-AT-KS', 'Trans-AT-KS',
                            'Trans-AT-KS', 'Trans-AT-KS', 'Trans-AT-KS']
        ks_subsubtypes_true = ['ARST', 'ALPHAME_EDB', 'ALPHAME', 'ALPHAME',
                               'EDB', 'BETA_MEDB', 'NON_ELONGATING_BETA_OH',
                               'EDB', 'NON_ELONGATING_BETA_L_OH', 'EDB',
                               'ALPHAME_BETAOH', 'BETA_L_OH', 'ALPHAME_EDB',
                               'NON_ELONGATING_DB']
        ks_subsubtypes = []
        for cds, cds_result in results_domains.cds_results.items():
            ks_subtypes += cds_result.ks_subtypes
            ks_subsubtypes += cds_result.ks_subsubtypes
        assert ks_subtypes == ks_subtypes_true
        assert ks_subsubtypes == ks_subsubtypes_true

    def test_kirromycin(self):
        filename = helpers.get_path_to_kirromycin_genbank()
        results = helpers.run_and_regenerate_results_for_module(filename, nrps_pks, self.options)
        assert len(results.domain_predictions) == 18
        domains = [('kirB_AMP-binding', 1), ('kirAII_PKS_KR', 1),
                   ('kirAII_PKS_KR', 2), ('kirAV_PKS_KR', 1),
                   ('kirAVI_PKS_AT', 2), ('kirAIV_PKS_KR', 3),
                   ('kirAI_PKS_KR', 1), ('kirAIV_PKS_KR', 2),
                   ('kirAIII_AMP-binding', 1), ('kirAVI_PKS_KR', 1),
                   ('kirAII_PKS_KR', 3), ('kirAIV_PKS_KR', 4),
                   ('kirAIV_PKS_KR', 1), ('kirCII_PKS_AT', 1),
                   ('kirCI_PKS_AT', 1), ('kirCI_PKS_AT', 2),
                   ('kirAV_PKS_KR', 2), ('kirAVI_PKS_AT', 1)]
        domain_names = ['nrpspksdomains_%s.%d' % dom for dom in domains]
        assert set(results.domain_predictions) == set(domain_names)
        assert set(results.domain_predictions[domain_names[0]]) == {"NRPSPredictor2"}
        pred = results.region_predictions[1][0]
        monomers = '(ccmal) + (Me-ccmal - ccmal - mal) + (gly - ccmal) + (ohmal - ohmal) + (ccmal) + (Me-ccmal - mal) + (gly)'
        assert pred.polymer == monomers
        results_domains = helpers.run_and_regenerate_results_for_module(filename, nrps_pks_domains, self.options)
        ks_subtypes = []
        ks_subtypes_true = ['Trans-AT-KS', 'Trans-AT-KS', 'Trans-AT-KS',
                            'Trans-AT-KS', 'Trans-AT-KS', 'Trans-AT-KS', 'Trans-AT-KS',
                            'Trans-AT-KS', 'Trans-AT-KS', 'Trans-AT-KS', 'Trans-AT-KS',
                            'Trans-AT-KS', 'Modular-KS', 'Modular-KS']

        ks_subsubtypes_true = ['ST', 'NON_ELONGATING_BETA_OH', 'DB', 'unknown variant',
                               'AA', 'unknown variant', 'TRANS_AT_PKS', 'TRANS_AT_PKS', 'unknown variant',
                               'BETA_OH_KETO', 'NON_ELONGATING_BETA_OH', 'DB']
        ks_subsubtypes = []
        for cds, cds_result in results_domains.cds_results.items():
            ks_subtypes += cds_result.ks_subtypes
            ks_subsubtypes += cds_result.ks_subsubtypes
        assert ks_subtypes == ks_subtypes_true
        assert ks_subsubtypes == ks_subsubtypes_true
