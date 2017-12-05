# License: GNU Affero General Public License v3 or later
# A copy of GNU AGPL v3 should have been included in this software package in LICENSE.txt.

# for test files, silence irrelevant and noisy pylint warnings
# pylint: disable=no-self-use,protected-access,missing-docstring

import unittest

from Bio.SeqFeature import FeatureLocation
from minimock import mock, restore

from antismash.common import subprocessing  # mocked, pylint: disable=unused-import
from antismash.common.test.helpers import DummyCDS, FakeHit
from antismash.modules.thiopeptides.specific_analysis import (
    Thiopeptide,
    predict_cleavage_site,
    result_vec_to_feature,
    acquire_rodeo_heuristics,
)


class TestThiopeptide(unittest.TestCase):
    def test_init(self):
        "Test Thiopeptide instantiation"
        thio = Thiopeptide(20, 31, 32, 51)
        assert isinstance(thio, Thiopeptide)
        self.assertEqual(20, thio.start)
        self.assertEqual(31, thio.end)
        self.assertEqual(32, thio.score)
        self.assertEqual('', thio.thio_type)
        self.assertEqual('', thio.core)
        self.assertEqual(51, thio.rodeo_score)
        with self.assertRaisesRegex(AssertionError, "calculating weight without a core"):
            print(thio.molecular_weight)

    def test_repr(self):
        "Test Thiopeptide representation"
        thio = Thiopeptide(20, 31, 32, 51)
        thio.core = "MAGIC"
        expected = "Thiopeptide(20..31, 32, 'MAGIC', '', -1(-1), , False, , )"
        self.assertEqual(expected, repr(thio))

    def test_core(self):
        "Test Thiopeptide.core"
        thio = Thiopeptide(20, 31, 32, 51)
        thio.core = "MAGICHAT"
        self.assertEqual('MAGICHAT', thio.core)
        assert thio.core_analysis
        thio.core = "MAGICXHAT"
        self.assertEqual('MAGICXHAT', thio.core)
        assert thio.core_analysis

    def test_macrocycle_size(self):
        "Test Thiopeptide._predict_macrocycle"
        thio = Thiopeptide(20, 31, 32, 51)
        thio.core = "SAAAAAAAASC"
        self.assertEqual('26-member', thio.macrocycle)
        thio.core = "SAAAAAAAAASSSSSS"
        self.assertEqual('35-member', thio.macrocycle)
        thio.core = "SDDDDDDDDDSC"
        self.assertEqual('29-member', thio.macrocycle)
        thio.core = "TAAASDDDDDDDDSCDD"
        self.assertEqual('26-member', thio.macrocycle)

    def test_weights(self):
        "Test Thiopeptide.alt_mature_weights"
        # includes all weights (even after maturation)
        thio = Thiopeptide(23, 42, 17, 51)
        thio.core = "MAGICCHATS"
        thio.amidation = True
        thio.thio_type = "Type-I"

        mat_weights = thio.mature_alt_weights
        weight = mat_weights[0]
        mono = mat_weights[1]
        alt = mat_weights[2:]
        self.assertAlmostEqual(weight, 1051.1623)
        self.assertAlmostEqual(mono, 1050.367788)
        self.assertAlmostEqual(alt, [1069.1823, 1087.2023])

        thio.core = "MAGICCHATS"
        thio.amidation = False
        thio.thio_type = "Type-I"
        # reset mature alt weights to calculate them again
        thio._mature_alt_weights = []

        mat_weights = thio.mature_alt_weights
        weight = mat_weights[0]
        mono = mat_weights[1]
        alt = mat_weights[2:]
        self.assertAlmostEqual(weight, 1121.1623)
        self.assertAlmostEqual(mono, 1120.367788)
        self.assertAlmostEqual(alt, [1139.1823, 1157.2023])

        thio.core = "MAGICCHATS"
        thio.amidation = True
        thio.thio_type = "Type-II"
        # reset mature alt weights to calculate them again
        thio._mature_alt_weights = []

        mat_weights = thio.mature_alt_weights
        weight = mat_weights[0]
        mono = mat_weights[1]
        alt = mat_weights[2:]
        self.assertAlmostEqual(weight, 1097.1623)
        self.assertAlmostEqual(mono, 1096.367788)
        self.assertAlmostEqual(alt, [1115.1823, 1133.2023])


class TestSpecificAnalysis(unittest.TestCase):
    def setUp(self):
        self.hmmpfam_return_vals = []
        mock('subprocessing.run_hmmpfam2', returns=self.hmmpfam_return_vals)

    def tearDown(self):
        restore()

    def test_predict_cleavage_site(self):
        "Test thiopeptides.predict_cleavage_site()"
        resvec = predict_cleavage_site('foo', 'bar', 51)
        self.assertEqual([None, None, None], resvec)
        fake_hit = FakeHit(24, 42, 17, 'fake')
        self.hmmpfam_return_vals.append([fake_hit])

        start, end, score = predict_cleavage_site('foo', 'bar', 15)

        self.assertEqual(24, start)
        self.assertEqual(28, end)
        self.assertEqual(17, score)

    def test_result_vec_to_feature(self):
        "Test thiopeptides.result_vec_to_features()"
        loc = FeatureLocation(0, 165, strand=1)
        orig_feature = DummyCDS(0, 165, locus_tag='FAKE0001')
        vec = Thiopeptide(17, 23, 42, 51)
        seq = 'SCTSSCTSS'
        vec.thio_type = 'Type-III'
        vec.core = seq
        vec.leader = "HEADHEADHEAD"
        orig_feature.translation = seq + vec.leader
        motif = result_vec_to_feature(orig_feature, vec)

        assert loc.start == motif.leader.start
        assert loc.start + (23 * 3) == motif.leader.end
        assert loc.strand == motif.leader.strand
        assert motif.type == 'CDS_motif'
        assert motif.peptide_type == "thiopeptide"
        assert orig_feature.locus_tag == motif.locus_tag
        assert motif.rodeo_score == 51
        assert motif.score == 42
        self.assertAlmostEqual(motif.molecular_weight, 861.9, places=1)
        assert motif.peptide_class == "Type III"

        assert motif.leader_seq == "HEADHEADHEAD"
        assert motif.leader.end == motif.core.start
        assert loc.end == motif.core.end
        assert loc.strand == motif.core.strand
        self.assertAlmostEqual(motif.monoisotopic_mass, 861.3, places=1)
        assert len(motif.alternative_weights) == 7
        for calc, expect in zip(motif.alternative_weights, [879.9, 897.9, 916.0,
                                                            934.0, 952.0, 970.0,
                                                            988.0]):
            self.assertAlmostEqual(calc, expect, places=1)
        assert not motif.amidation
        assert not motif.macrocycle
        assert not motif.cleaved_residues
        assert motif.core_features == "Central ring: pyridine trisubstituted"
        assert motif.core_seq == "SCTSSCTSS"

    def test_acquire_rodeo_heuristics(self):
        """Test thiopeptides.acquire_rodeo_heuristics()"""
        leader = "MSDITASRVESLDLQDLDLSELTVTSLRDTVALPENGA"
        core = "SWGSCSCQASSSCAQPQDM"
        domains = [
            'LANC_like',
            'Pkinase',
            'Lant_dehyd_N',
            'Lant_dehyd_C',
            'Lant_dehydr_C',
            'YcaO',
            'PF00881',
            'Lant_dehydr_C'
        ]
        expected_score = 15
        expected_tabs = [1, 1, 1, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 1, 0, 1, 1, 0, 0, 0, 0, 1]
        print(expected_tabs)
        score, tabs = acquire_rodeo_heuristics(leader, core, domains)
        print(tabs)
        self.assertEqual(expected_score, score)
        self.assertEqual(expected_tabs, tabs)

        core = "SWGSCSCQASSSCAQPQDMX"
        score, tabs = acquire_rodeo_heuristics(leader, core, domains)
        self.assertEqual(expected_score, score)
        self.assertEqual(expected_tabs, tabs)
