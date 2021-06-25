# License: GNU Affero General Public License v3 or later
# A copy of GNU AGPL v3 should have been included in this software package in LICENSE.txt.

# for test files, silence irrelevant and noisy pylint warnings
# pylint: disable=no-self-use,protected-access,missing-docstring

import unittest
from unittest.mock import patch

from Bio.SeqUtils.ProtParam import ProteinAnalysis

from antismash.common import subprocessing
from antismash.common.test.helpers import DummyCDS, FakeHSP

from antismash.modules.lassopeptides.specific_analysis import (
    Lassopeptide,
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
