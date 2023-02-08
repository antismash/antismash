# License: GNU Affero General Public License v3 or later
# A copy of GNU AGPL v3 should have been included in this software package in LICENSE.txt.

# for test files, silence irrelevant and noisy pylint warnings
# pylint: disable=use-implicit-booleaness-not-comparison,protected-access,missing-docstring

import unittest
from unittest.mock import patch

from Bio.SeqUtils.ProtParam import ProteinAnalysis

from antismash.common import subprocessing
from antismash.common.test.helpers import DummyCDS, FakeHit, DummyRecord
from antismash.modules.lanthipeptides import specific_analysis as lanthi
from antismash.modules.lanthipeptides.specific_analysis import (
    Lanthipeptide,
    LanthiResults,
    predict_cleavage_site,
    result_vec_to_feature,
    CleavageSiteHit,
    run_lanthipred
)


class TestLanthipeptide(unittest.TestCase):
    def setUp(self):
        self.lant = Lanthipeptide(CleavageSiteHit(42, 17, 'Class-I'), 23, "HEAD", "MAGICHAT")

    def test_init(self):
        "Test Lanthipeptide instantiation"
        assert isinstance(self.lant, Lanthipeptide)
        assert self.lant.end == 42
        assert self.lant.score == 17
        assert self.lant.lantype == "Class-I"
        assert self.lant.leader == "HEAD"
        assert self.lant.core == "MAGICHAT"
        assert self.lant.core_analysis

    def test_repr(self):
        "Test Lanthipeptide representation"
        expected = "Lanthipeptide(..42, 17.0, 'Class-I', 'MAGICHAT', -1, "  # skip weights
        assert repr(self.lant).startswith(expected)

    def test_core_ignore_invalid(self):
        "Test Lanthipeptide.core ignores invalid amino acids"
        with self.assertRaisesRegex(AssertionError, "calculating weight without a core"):
            self.lant.core = ""

        self.lant.core = "MAGICXHAT"
        assert self.lant.core == "MAGICXHAT"
        assert self.lant.core_analysis

    def test_number_of_lan_bridges(self):
        "Test Lanthipeptide.number_of_lan_bridges"
        self.lant.core = "MAGICHAT"
        assert self.lant.number_of_lan_bridges == 1
        self.lant.core = "MAGICHATCS"
        assert self.lant.number_of_lan_bridges == 2
        self.lant.core = "MAGICHATCSS"
        assert self.lant.number_of_lan_bridges == 2
        self.lant.core = "MAGICHATCT"
        assert self.lant.number_of_lan_bridges == 2
        self.lant.core = "MAGICHATCCS"
        assert self.lant.number_of_lan_bridges == 2

    def test_monoisotopic_mass(self):
        "Test Lanthipeptide.monoisotopic_mass"
        assert self.lant.core == "MAGICHAT"
        analysis = ProteinAnalysis("MAGICHAT", monoisotopic=True)
        weight = analysis.molecular_weight()
        # Thr is assumed to be dehydrated
        weight -= 18
        self.assertAlmostEqual(weight, self.lant.monoisotopic_mass)
        self.assertAlmostEqual(weight, self.lant._monoisotopic_weight)

    def test_molecular_weight(self):
        "Test Lanthipeptide.molecular_weight"
        assert self.lant.core == "MAGICHAT"
        analysis = ProteinAnalysis("MAGICHAT", monoisotopic=False)
        weight = analysis.molecular_weight()
        # Thr is assumed to be dehydrated
        weight -= 18.02
        self.assertAlmostEqual(weight, self.lant.molecular_weight)
        self.assertAlmostEqual(weight, self.lant._weight)

    def test_alternative_weights(self):
        "Test Lanthipeptide.alt_weights"
        self.lant.core = "MAGICHATS"
        analysis = ProteinAnalysis("MAGICHATS", monoisotopic=False)
        weight = analysis.molecular_weight()
        # One Ser/Thr is assumed to be dehydrated, but not the other
        weight -= 18.02
        self.assertEqual([weight], self.lant.alternative_weights)


class TestSpecificAnalysis(unittest.TestCase):
    def test_predict_cleavage_site(self):
        with patch.object(subprocessing, "run_hmmpfam2", return_value=[]):
            assert predict_cleavage_site('foo', 'bar') is None

        fake_hit = FakeHit(24, 42, 17, 'fake')
        with patch.object(subprocessing, "run_hmmpfam2", return_value=[[fake_hit]]):
            res = predict_cleavage_site('foo', 'bar')
        self.assertEqual(42, res.end)
        self.assertEqual(17, res.score)
        self.assertEqual('fake', res.lantype)

    def test_result_vec_to_features(self):
        "Test lanthipeptides.result_vec_to_features()"
        orig_feature = DummyCDS(0, 165)
        orig_feature.locus_tag = 'FAKE0001'
        seq = "TAILTAILTAILTAILTAILTAILTAILTAILTAILCC"
        vec = Lanthipeptide(CleavageSiteHit(23, 42, 'Class-I'), 23, "HEADHEADHEAD", seq)
        motif = result_vec_to_feature(orig_feature, vec)
        assert motif.location.start == 0
        assert motif.location.end == 165
        assert motif.location.strand == 1
        assert motif.locus_tag == f"{orig_feature.locus_tag}_{motif.peptide_class}"
        assert motif.score == 42
        assert motif.detailed_information.rodeo_score == 23
        assert motif.peptide_class == "lanthipeptide"
        assert motif.peptide_subclass == "Class I"
        assert motif.detailed_information.lan_bridges == 2
        self.assertAlmostEqual(motif.molecular_weight, 3648.6, places=1)
        self.assertAlmostEqual(motif.monoisotopic_mass, 3646.3, places=1)
        for calc, expected in zip(motif.alternative_weights,
                                  [3666.6, 3684.6, 3702.7, 3720.7, 3738.7, 3756.7, 3774.7]):
            self.assertAlmostEqual(calc, expected, places=1)
        assert motif.core == seq
        assert motif.leader == vec.leader


@patch.object(lanthi, "run_rodeo", return_value=20)
class TestNoCores(unittest.TestCase):
    """ Ensure that cleavage sites that result in no prepeptide core don't
        generate hits """
    def setUp(self):
        self.domains = ['Condensation', 'AMP-binding', 'PP-binding', 'Condensation',
                        'A-OX', 'PP-binding', 'PKS_AT', 'hyb_KS', 'PP-binding',
                        'PF00561', 'Peptidase_S9', 'LANC_like', 'Pkinase', 'p450',
                        'Abi', 'DUF4135', 'Condensation', 'AMP-binding', 'PP-binding',
                        'Condensation', 'AMP-binding', 'PP-binding', 'PF07366',
                        'adh_short', 'adh_short_C2', 'p450', 'PF12697', 'PP-binding',
                        'PKS_AT', 'mod_KS', 'hyb_KS', 'adh_short', 'PP-binding',
                        'PKS_AT', 'mod_KS', 'adh_short', 'PP-binding', 'PKS_AT',
                        'mod_KS', 'adh_short', 'NAD_binding_4', 'Condensation',
                        'PP-binding']
        self.cds = DummyCDS(38, 48)
        # an all_orf detected gene at 7922405:7922549 in CP013129.1, with a C appended
        self.cds.translation = "VCGPRDHGRQRTSAPHAFHSDSMDRAASRPVEYGDYSGSPLSQGLGGC"

    def test_prediction_with_no_core(self, _patched_rodeo):
        # the real cleavage site result (+1 at end for the C)
        cleavage_result = CleavageSiteHit(end=48, score=-6.8, lantype="Class-II")
        with patch.object(lanthi, "predict_cleavage_site", return_value=cleavage_result):
            for part in ["I", "II"]:
                assert run_lanthipred(None, self.cds, f"Class-{part}", self.domains) is None

    def test_prediction_with_core_class1(self, _patched_rodeo):
        # the cleavage result adjusted to leave at least one amino in core
        cleavage_result = CleavageSiteHit(end=40, score=-6.8, lantype="Class-I")
        with patch.object(lanthi, "predict_cleavage_site", return_value=cleavage_result):
            results = run_lanthipred(DummyRecord(features=[self.cds]),
                                     self.cds, "Class-I", self.domains)
        assert results
        assert str(results).startswith("Lanthipeptide(..40, -6.8, 'Class-I', 'LSQGLGGC', 1, 715")

    def test_prediction_with_core_class2(self, _patched_rodeo):
        # the cleavage result adjusted to leave at least one amino in core
        cleavage_result = CleavageSiteHit(end=40, score=-6.8, lantype="Class-II")
        with patch.object(lanthi, "predict_cleavage_site", return_value=cleavage_result):
            results = run_lanthipred(DummyRecord(features=[self.cds]),
                                     self.cds, "Class-II", self.domains)
        assert results is not None
        assert str(results).startswith("Lanthipeptide(..40, -6.8, 'Class-II', 'LSQGLGGC', 1, 715")


class TestCDSDuplication(unittest.TestCase):
    def setUp(self):
        self.record = DummyRecord(seq="A" * 100)
        self.cds = DummyCDS()

    def test_existing(self):
        self.record.add_cds_feature(self.cds)
        results = LanthiResults(self.record.id)
        results.add_cds(self.cds)
        assert len(self.record.get_cds_features()) == 1
        results.add_to_record(self.record)
        assert len(self.record.get_cds_features()) == 1

    def test_new(self):
        results = LanthiResults(self.record.id)
        results.add_cds(self.cds)
        assert len(self.record.get_cds_features()) == 0
        results.add_to_record(self.record)
        assert len(self.record.get_cds_features()) == 1

    def test_double_protocluster(self):
        results = LanthiResults(self.record.id)
        assert len(results._new_cds_features) == 0
        results.add_cds(self.cds)
        assert len(results._new_cds_features) == 1
        results.add_cds(self.cds)
        assert len(results._new_cds_features) == 1
        results.add_cds(DummyCDS(locus_tag="different"))
        assert len(results._new_cds_features) == 2
