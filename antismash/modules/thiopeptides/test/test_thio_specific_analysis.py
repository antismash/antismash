# License: GNU Affero General Public License v3 or later
# A copy of GNU AGPL v3 should have been included in this software package in LICENSE.txt.

# for test files, silence irrelevant and noisy pylint warnings
# pylint: disable=use-implicit-booleaness-not-comparison,protected-access,missing-docstring,unbalanced-tuple-unpacking

import unittest
from unittest.mock import patch

from Bio.SeqFeature import FeatureLocation

from antismash.common import all_orfs, subprocessing
from antismash.common.test.helpers import (
    DummyCDS,
    DummyProtocluster,
    DummyRecord,
    FakeHit,
)
from antismash.modules.thiopeptides.specific_analysis import (
    Thiopeptide,
    ThioResults,
    find_unannotated_candidates,
    predict_cleavage_site,
    result_vec_to_feature,
)


class TestThiopeptide(unittest.TestCase):
    def test_init(self):
        "Test Thiopeptide instantiation"
        thio = Thiopeptide(31, 32, 51)
        assert isinstance(thio, Thiopeptide)
        self.assertEqual(31, thio.end)
        self.assertEqual(32, thio.score)
        self.assertEqual('', thio.thio_type)
        self.assertEqual('', thio.core)
        self.assertEqual(51, thio.rodeo_score)
        with self.assertRaisesRegex(AssertionError, "calculating weight without a core"):
            print(thio.molecular_weight)

    def test_repr(self):
        "Test Thiopeptide representation"
        thio = Thiopeptide(31, 32, 51)
        thio.core = "MAGIC"
        thio.thio_type = "Type I"
        thio.c_cut = "DUMMYC"
        thio.macrocycle = "dummy macro"
        expected = ("Thiopeptide(..31, 32, 'MAGIC', 'Type I', -1.0(-1.0), dummy macro, False, "
                    "Central ring: pyridine tetrasubstituted (hydroxyl group present); second macrocycle, "
                    "DUMMYC)")
        self.assertEqual(expected, repr(thio))

    def test_core(self):
        "Test Thiopeptide.core"
        thio = Thiopeptide(31, 32, 51)
        thio.core = "MAGICHAT"
        self.assertEqual('MAGICHAT', thio.core)
        assert thio.core_analysis
        thio.core = "MAGICXHAT"
        self.assertEqual('MAGICXHAT', thio.core)
        assert thio.core_analysis

    def test_macrocycle_size(self):
        "Test Thiopeptide._predict_macrocycle"
        thio = Thiopeptide(31, 32, 51)
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
        thio = Thiopeptide(42, 17, 51)
        thio.core = "MAGICCHATS"
        thio.amidation = True
        thio.thio_type = "Type I"

        mat_weights = thio.mature_alt_weights
        weight = mat_weights[0]
        mono = mat_weights[1]
        alt = mat_weights[2:]
        self.assertAlmostEqual(weight, 1051.1623)
        self.assertAlmostEqual(mono, 1050.367788)
        self.assertAlmostEqual(alt, [1069.1823, 1087.2023])

        thio.core = "MAGICCHATS"
        thio.amidation = False
        thio.thio_type = "Type I"
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
        thio.thio_type = "Type II"
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
    @patch.object(subprocessing, 'run_hmmpfam2', return_value=[])
    def test_predict_cleavage_site(self, mocked_run):
        "Test thiopeptides.predict_cleavage_site()"
        resvec = predict_cleavage_site('foo', 'bar', 51)
        assert resvec == (None, 0.)
        fake_hit = FakeHit(24, 42, 17, 'fake')
        mocked_run.return_value = [[fake_hit]]

        end, score = predict_cleavage_site('foo', 'bar', 15)

        self.assertEqual(28, end)
        self.assertEqual(17, score)

    def test_result_vec_to_feature(self):
        "Test thiopeptides.result_vec_to_features()"
        loc = FeatureLocation(0, 66, strand=1)
        orig_feature = DummyCDS(0, 66, locus_tag='FAKE0001')
        vec = Thiopeptide(23, 42, 51)
        seq = 'SCTSSCTSS'
        vec.thio_type = 'Type III'
        vec.core = seq
        vec.leader = "HEADHEADHEAD"
        orig_feature.translation = seq + vec.leader
        motif = result_vec_to_feature(orig_feature, vec)

        leader, core = motif.to_biopython()

        assert loc.start == leader.location.start
        assert loc.start + (12 * 3) == leader.location.end
        assert loc.strand == leader.location.strand
        assert motif.type == 'CDS_motif'
        assert motif.peptide_class == "thiopeptide"
        assert motif.peptide_subclass == "Type III"
        assert motif.locus_tag == f"{orig_feature.locus_tag}_{motif.peptide_class}"
        assert motif.detailed_information.rodeo_score == 51
        assert motif.score == 42
        self.assertAlmostEqual(motif.molecular_weight, 861.9, places=1)

        assert motif.leader == "HEADHEADHEAD"
        assert leader.location.end == core.location.start
        assert loc.end == core.location.end
        assert loc.strand == core.location.strand
        self.assertAlmostEqual(motif.monoisotopic_mass, 861.3, places=1)
        assert len(motif.alternative_weights) == 7
        for calc, expect in zip(motif.alternative_weights, [879.9, 897.9, 916.0,
                                                            934.0, 952.0, 970.0,
                                                            988.0]):
            self.assertAlmostEqual(calc, expect, places=1)
        assert not motif.detailed_information.amidation
        assert not motif.detailed_information.macrocycle
        assert not motif.tail
        assert motif.detailed_information.core_features == "Central ring: pyridine trisubstituted"
        assert motif.core == "SCTSSCTSS"

    @patch.object(all_orfs, "find_all_orfs", return_value=[])
    @patch.object(subprocessing, "run_hmmscan", return_value=[])
    def test_find_unannotated(self, mocked_hmmscan, mocked_all_orfs):  # inverse order of decorators
        proto = DummyProtocluster(product="thiopeptide")
        record = DummyRecord()
        assert not find_unannotated_candidates(record, proto)
        assert mocked_all_orfs.called
        assert not mocked_hmmscan.called


class TestCDSDuplication(unittest.TestCase):
    def setUp(self):
        self.record = DummyRecord(seq="A" * 100)
        self.cds = DummyCDS()

    def test_existing(self):
        self.record.add_cds_feature(self.cds)
        results = ThioResults(self.record.id)
        results.add_cds(1, self.cds)
        assert len(self.record.get_cds_features()) == 1
        results.add_to_record(self.record)
        assert len(self.record.get_cds_features()) == 1

    def test_new(self):
        results = ThioResults(self.record.id)
        results.add_cds(1, self.cds)
        assert len(self.record.get_cds_features()) == 0
        results.add_to_record(self.record)
        assert len(self.record.get_cds_features()) == 1

    def test_double_protocluster(self):
        results = ThioResults(self.record.id)
        assert len(results._cds_features[1]) == 0
        results.add_cds(1, self.cds)
        assert len(results._cds_features[1]) == 1
        results.add_cds(1, self.cds)
        assert len(results._cds_features[1]) == 1
        assert len(results._cds_features) == 1
        results.add_cds(2, DummyCDS(start=50, end=80, locus_tag="different"))
        assert len(results._cds_features) == 2
