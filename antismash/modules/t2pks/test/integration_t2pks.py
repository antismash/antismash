# License: GNU Affero General Public License v3 or later
# A copy of GNU AGPL v3 should have been included in this software package in LICENSE.txt.

# for test files, silence irrelevant and noisy pylint warnings
# pylint: disable=use-implicit-booleaness-not-comparison,protected-access,missing-docstring

import unittest

import antismash
from antismash.common import path
from antismash.common.test import helpers
from antismash.config import build_config, destroy_config
from antismash.modules import t2pks


class T2PKSTest(unittest.TestCase):
    def setUp(self):
        self.options = build_config(["--minimal", "--enable-t2pks", "--enable-html"], isolated=True,
                                    modules=antismash.get_all_modules())

    def tearDown(self):
        destroy_config()

    def test_reuse(self):
        test_file = path.get_full_path(__file__, 'data', 'NC_003888.3.cluster011.gbk')
        results = helpers.run_and_regenerate_results_for_module(test_file, t2pks, self.options)
        assert list(results.cluster_predictions) == [1]
        pred = results.cluster_predictions[1]
        assert pred.starter_units == [t2pks.results.Prediction('acetyl-CoA', 0., 0.)]
        assert pred.malonyl_elongations == [t2pks.results.Prediction('7', 743.5, 1.2e-226)]
        assert pred.product_classes == {'benzoisochromanequinone'}
        assert list(pred.molecular_weights) == ['acetyl-CoA_7']
        self.assertAlmostEqual(pred.molecular_weights['acetyl-CoA_7'], 342.3845)

    def test_full_blastp_use(self):
        test_file = path.get_full_path(__file__, 'data', 'GQ409537.1.gbk')
        results = helpers.run_and_regenerate_results_for_module(test_file, t2pks, self.options)
        assert list(results.cluster_predictions) == [1]
        pred = results.cluster_predictions[1]
        assert pred.starter_units == [t2pks.results.Prediction('malonamyl-CoA', 2319., 0.)]
        assert pred.malonyl_elongations == [t2pks.results.Prediction('8|9', 661.0, 1.3e-201)]
        assert pred.product_classes == {'angucycline', 'tetracycline', 'aureolic acid',
                                        'anthracycline', 'tetracenomycin'}
        assert set(pred.molecular_weights) == {'malonamyl-CoA_8', 'malonamyl-CoA_9'}
        self.assertAlmostEqual(pred.molecular_weights['malonamyl-CoA_8'], 654.63434)
        self.assertAlmostEqual(pred.molecular_weights['malonamyl-CoA_9'], 696.67102)
