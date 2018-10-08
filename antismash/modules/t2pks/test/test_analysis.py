# License: GNU Affero General Public License v3 or later
# A copy of GNU AGPL v3 should have been included in this software package in LICENSE.txt.

# for test files, silence irrelevant and noisy pylint warnings
# pylint: disable=no-self-use,protected-access,missing-docstring

import unittest

from minimock import mock, restore

from antismash.common import path, record_processing
from antismash.common.hmmscan_refinement import HMMResult
from antismash.common.secmet.features.cluster import Cluster, FeatureLocation
from antismash.modules.t2pks import t2pks_analysis
from antismash.modules.t2pks.results import ClusterPrediction, Prediction


class TestCoelicolorAnalysis(unittest.TestCase):
    def setUp(self):
        test_file = path.get_full_path(__file__, 'data', 'NC_003888.3.cluster011.gbk')
        self.record = record_processing.parse_input_sequence(test_file)[0]
        self.cluster = Cluster(FeatureLocation(0, len(self.record.seq)),
                               surrounding_location=FeatureLocation(0, len(self.record.seq)),
                               cutoff=20, neighbourhood_range=0, tool="test", product="t2pks",
                               detection_rule="dummy rule")
        self.record.add_cluster(self.cluster)
        self.record.create_superclusters()
        self.record.create_regions()
        hmm_results = {'SCO5072': [HMMResult("KR", 1, 265, evalue=3.1e-49, bitscore=159.4)],
                       'SCO5079': [HMMResult("DIMER", 4, 293, evalue=8.7e-131, bitscore=426.8)],
                       'SCO5080': [HMMResult("OXY", 8, 377, evalue=2.1e-14, bitscore=44.7)],
                       'SCO5086': [HMMResult("KR_C9", 0, 261, evalue=1.9e-134, bitscore=438.4)],
                       'SCO5087': [HMMResult("KS", 44, 463, evalue=3.5e-234, bitscore=768.6)],
                       'SCO5088': [HMMResult("CLF_7", 1, 401, evalue=1.2e-226, bitscore=743.5)],
                       'SCO5089': [HMMResult("ACP", 4, 86, evalue=5e-36, bitscore=114.2)],
                       'SCO5090': [HMMResult("CYC_C7-C12", 1, 312, evalue=7.8e-124, bitscore=404)],
                       'SCO5091': [HMMResult("CYC_C5-C14", 3, 297, evalue=4.4e-143, bitscore=467.3)],
                       'SCO5094': [HMMResult("MET", 40, 155, evalue=9.8e-11, bitscore=32.7)],
                       'SCO5097': [HMMResult("KR", 3, 247, evalue=3.3e-40, bitscore=129.8)],
                       }
        mock("t2pks_analysis.run_t2pks_hmmscan", returns=hmm_results)
        mock("t2pks_analysis.run_starter_unit_blastp", returns={})

    def tearDown(self):
        restore()

    def test_coelicolor_c11(self):
        results = t2pks_analysis.analyse_cluster(self.cluster, self.record)
        assert isinstance(results, ClusterPrediction)
        assert results.product_classes == {"benzoisochromanequinone"}
        assert results.starter_units == [Prediction("acetyl", 0., 0.)]
        assert results.malonyl_elongations == [Prediction("7", 743.5, 1.2e-226)]
        assert list(results.molecular_weights) == ["acetyl_7"]
        self.assertAlmostEqual(results.molecular_weights["acetyl_7"], 451.22572)
