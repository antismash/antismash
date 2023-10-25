# License: GNU Affero General Public License v3 or later
# A copy of GNU AGPL v3 should have been included in this software package in LICENSE.txt.

# for test files, silence irrelevant and noisy pylint warnings
# pylint: disable=use-implicit-booleaness-not-comparison,protected-access,missing-docstring

import unittest
from unittest.mock import patch

from Bio.SeqUtils.ProtParam import ProteinAnalysis

from antismash.common import subprocessing
from antismash.common.secmet.test.helpers import DummyProtocluster
from antismash.common.secmet.features import Prepeptide
from antismash.common.secmet.locations import FeatureLocation
from antismash.common.test.helpers import DummyCDS, DummyRecord, FakeHSP

from antismash.modules.lassopeptides import specific_analysis as module
from antismash.modules.lassopeptides.specific_analysis import (
    Lassopeptide,
    LassoResults,
    predict_cleavage_site,
    result_vec_to_motif,
)


class TestLassopeptide(unittest.TestCase):
    def test_init(self):
        "Test Lassopeptide instantiation"
        lasso = Lassopeptide(42, 17, 51, "", "")
        self.assertTrue(isinstance(lasso, Lassopeptide))
        self.assertEqual(42, lasso.end)
        self.assertEqual(17, lasso.score)
        self.assertEqual('Class II', lasso._lassotype)
        self.assertEqual('', lasso.core)
        self.assertEqual(51, lasso.rodeo_score)
        with self.assertRaises(ValueError):
            print(lasso.molecular_weight)

    def test_repr(self):
        "Test Lassopeptide representation"
        lasso = Lassopeptide(42, 17, 51, "lead", "ABCDEFGHIJKLM")
        expected = "Lassopeptide(..42, 17, 'Class II', 'ABCDEFGHIJKLM', 0, -1(-1), , )"
        self.assertEqual(expected, repr(lasso))

    def test_core(self):
        "Test Lassopeptide.core"
        lasso = Lassopeptide(42, 17, 51, "lead", "MAGICHAT")
        assert lasso.core == "MAGICHAT"
        self.assertAlmostEqual(lasso.monoisotopic_mass, 784.3, delta=1)

    def test_core_ignore_invalid(self):
        "Test Lassopeptide.core ignores invalid amino acids"
        lasso = Lassopeptide(42, 17, 51, "lead", "MAGICXHAT")
        assert lasso.core == "MAGICXHAT"
        self.assertAlmostEqual(lasso.monoisotopic_mass, 784.3, delta=1)

    def test_number_bridges_class(self):
        "Test Lassopeptide.number_bridges"
        lasso = Lassopeptide(42, 20, 51, "lead", "MAGIXHAT")
        self.assertEqual(0, lasso.number_bridges)
        self.assertEqual('Class II', lasso.lasso_class)
        lasso.core = "CAGIMHAC"
        self.assertEqual(1, lasso.number_bridges)
        self.assertEqual('Class III', lasso.lasso_class)
        lasso.core = "CCCC"
        self.assertEqual(2, lasso.number_bridges)
        self.assertEqual('Class I', lasso.lasso_class)

    def test_monoisotopic_mass(self):
        "Test Lassopeptide.monoisotopic_mass"
        lasso = Lassopeptide(42, 17, 51, "lead", "MAGICHATTIP")
        lasso.c_cut = "TIP"
        analysis = ProteinAnalysis("MAGICHATTIP", monoisotopic=True)
        cut_analysis = ProteinAnalysis("MAGICHAT", monoisotopic=True)
        mw = analysis.molecular_weight() - 18.02
        cut_mw = cut_analysis.molecular_weight() - 18.02
        self.assertAlmostEqual(mw, lasso.monoisotopic_mass)
        self.assertAlmostEqual(cut_mw, lasso.cut_mass)

    def test_molecular_weight(self):
        "Test Lassopeptide.molecular_weight"
        lasso = Lassopeptide(42, 17, 51, "lead", "MAGICHATTIP")
        lasso.c_cut = "TIP"
        analysis = ProteinAnalysis("MAGICHATTIP", monoisotopic=False)
        cut_analysis = ProteinAnalysis("MAGICHAT", monoisotopic=False)
        mw = analysis.molecular_weight() - 18.02
        cut_mw = cut_analysis.molecular_weight() - 18.02
        self.assertAlmostEqual(mw, lasso.molecular_weight)
        self.assertAlmostEqual(cut_mw, lasso.cut_weight)

    def test_macrolactam(self):
        "Test Lassopeptide.macrolactam"
        lasso = Lassopeptide(42, 20, 51, "lead", "GAAAAADLLLLLLLLL")
        self.assertEqual("GAAAAAD", lasso.macrolactam)


class TestSpecificAnalysis(unittest.TestCase):
    class FakeHit:  # pylint: disable=too-few-public-methods
        def __init__(self, start, end, score, desc):
            self.hsps = [FakeHSP(start, end, score)]
            self.description = desc

        def __iter__(self):
            return iter(self.hsps)

    def test_predict_cleavage_site(self):
        "Test lassopeptides.predict_cleavage_site()"
        with patch.object(subprocessing, "run_hmmpfam2", return_value=[]):
            assert predict_cleavage_site('foo', 'bar', 51) == (None, 0.)

        fake_hit = self.FakeHit(24, 42, 17, 'Class II')

        with patch.object(subprocessing, "run_hmmpfam2", return_value=[[fake_hit, fake_hit]]):
            assert predict_cleavage_site('foo', 'bar', 51) == (None, 0.)
            assert predict_cleavage_site('foo', 'bar', 15) == (41, 17)

    def test_result_vec_to_features(self):
        "Test lassopeptides.result_vec_to_features()"
        orig_feature = DummyCDS(0, 165, strand=1, locus_tag='FAKE0001')
        seq = "TAILTAILTAILTAILTAILTAILTAILTAILTAILCCTIP"
        vec = Lassopeptide(23, 42, 51, "HEADHEADHEAD", seq)
        vec.c_cut = "TIP"
        motif = result_vec_to_motif(orig_feature, vec)
        assert motif.locus_tag == f"{orig_feature.locus_tag}_{motif.peptide_class}"
        assert motif.location == orig_feature.location
        assert motif.leader == "HEADHEADHEAD"
        assert motif.tail == "TIP"
        assert motif.core == seq[:-len(motif.tail)]
        assert motif.detailed_information.num_bridges == 1
        assert motif.peptide_subclass == "Class III"
        assert motif.detailed_information.rodeo_score == 51
        assert motif.score == 42
        assert motif.detailed_information.macrolactam == ""
        self.assertAlmostEqual(motif.monoisotopic_mass, 4103.5, delta=1)
        self.assertAlmostEqual(motif.molecular_weight, 4106.1, delta=1)
        self.assertAlmostEqual(motif.detailed_information.cut_mass, 3792.3, delta=1)
        self.assertAlmostEqual(motif.detailed_information.cut_weight, 3794.7, delta=1)


class TestCDSDuplication(unittest.TestCase):
    def setUp(self):
        self.record = DummyRecord(seq="A" * 100)
        self.cds = DummyCDS()

    def test_existing(self):
        self.record.add_cds_feature(self.cds)
        results = LassoResults(self.record.id)
        results.add_cds(self.cds)
        assert len(self.record.get_cds_features()) == 1
        results.add_to_record(self.record)
        assert len(self.record.get_cds_features()) == 1

    def test_new(self):
        results = LassoResults(self.record.id)
        results.add_cds(self.cds)
        assert len(self.record.get_cds_features()) == 0
        results.add_to_record(self.record)
        assert len(self.record.get_cds_features()) == 1

    def test_double_protocluster(self):
        results = LassoResults(self.record.id)
        assert len(results._new_cds_features) == 0
        results.add_cds(self.cds)
        assert len(results._new_cds_features) == 1
        results.add_cds(self.cds)
        assert len(results._new_cds_features) == 1
        results.add_cds(DummyCDS(locus_tag="different"))
        assert len(results._new_cds_features) == 2


class TestMotifDuplication(unittest.TestCase):
    def test_shared_neighbourhood(self):
        record = DummyRecord(seq="A" * 120)
        overlapping = [
            DummyProtocluster(0, 50, core_start=10, core_end=20, product="lassopeptide"),
            DummyProtocluster(30, 80, core_start=60, core_end=70, product="lassopeptide"),
        ]
        for overlap in overlapping:
            record.add_protocluster(overlap)
        standalone = DummyProtocluster(100, 120, core_start=105, core_end=110, product="lassopeptide")
        record.add_protocluster(standalone)

        record.add_cds_feature(DummyCDS(start=30, end=36, locus_tag="shared"))
        record.add_cds_feature(DummyCDS(start=105, end=110, locus_tag="other"))
        # all the protoclusters need to have those cdses, otherwise there'll be problems later
        assert len(record.get_protoclusters()) == 3
        for proto in record.get_protoclusters():
            assert len(proto.cds_children) == 1, proto

        prepeptide = Prepeptide(FeatureLocation(55, 58, 1), "lassopeptide", "CORE", "shared", "lassopeptides")
        other = Prepeptide(FeatureLocation(155, 158, 1), "lassopeptide", "CORE",  "other", "lassopeptides")

        with patch.object(module, "run_lassopred", side_effect=[prepeptide, prepeptide, other]):
            with patch.object(module.all_orfs, "find_all_orfs", return_value=[]):
                results = module.specific_analysis(record)
        # both of the overlapping protoclusters should have references to it
        for cluster in overlapping:
            assert results.clusters[cluster.get_protocluster_number()] == {prepeptide.get_name()}
        # and the non-overlapping should not
        assert results.clusters[standalone.get_protocluster_number()] == {other.get_name()}
        # but the motif should not duplicate by name
        assert set(results.motifs_by_locus) == {prepeptide.get_name(), other.get_name()}
        assert len(results.motifs_by_locus[prepeptide.get_name()]) == 1
