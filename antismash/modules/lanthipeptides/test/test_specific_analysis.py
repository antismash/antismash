# License: GNU Affero General Public License v3 or later
# A copy of GNU AGPL v3 should have been included in this software package in LICENSE.txt.

import unittest

from Bio.SeqUtils.ProtParam import ProteinAnalysis
from minimock import mock, restore, TraceTracker

# pylint: disable=unused-import
from antismash.common import subprocessing # used in mocks
# pylint: enable=unused-import
from antismash.common.test.helpers import DummyCDS
from antismash.modules.lanthipeptides.specific_analysis import (
    Lanthipeptide,
    predict_cleavage_site,
    result_vec_to_feature,
)

class TestLanthipeptide(unittest.TestCase):
    def test_init(self):
        "Test Lanthipeptide instantiation"
        lant = Lanthipeptide(23, 42, 17, 23, 'Class-I')
        self.assertTrue(isinstance(lant, Lanthipeptide))
        self.assertEqual(23, lant.start)
        self.assertEqual(42, lant.end)
        self.assertEqual(17, lant.score)
        self.assertEqual('Class-I', lant.lantype)
        self.assertEqual('', lant.core)
        with self.assertRaises(ValueError):
            dummy = lant.molecular_weight

    def test_repr(self):
        "Test Lanthipeptide representation"
        lant = Lanthipeptide(23, 42, 17, 23, 'Class-I')
        expected = "Lanthipeptide(23..42, 17, 'Class-I', '', -1, -1(-1))"
        self.assertEqual(expected, repr(lant))

    def test_core(self):
        "Test Lanthipeptide.core"
        lant = Lanthipeptide(23, 42, 17, 23, 'Class-I')
        self.assertEqual('', lant.core)
        self.assertFalse(hasattr(lant, 'core_analysis'))
        lant.core = "MAGICHAT"
        self.assertEqual('MAGICHAT', lant.core)
        self.assertTrue(hasattr(lant, 'core_analysis'))

    def test_core_ignore_invalid(self):
        "Test Lanthipeptide.core ignores invalid amino acids"
        lant = Lanthipeptide(23, 42, 17, 23, 'Class-I')
        self.assertEqual('', lant.core)
        self.assertFalse(hasattr(lant, 'core_analysis'))
        lant.core = "MAGICXHAT"
        self.assertEqual('MAGICXHAT', lant.core)
        self.assertTrue(hasattr(lant, 'core_analysis'))

    def test_number_of_lan_bridges(self):
        "Test Lanthipeptide.number_of_lan_bridges"
        lant = Lanthipeptide(23, 42, 17, 23, 'Class-I')
        lant.core = "MAGICHAT"
        self.assertEqual(1, lant.number_of_lan_bridges)
        lant.core = "MAGICHATCS"
        self.assertEqual(2, lant.number_of_lan_bridges)
        lant.core = "MAGICHATCSS"
        self.assertEqual(2, lant.number_of_lan_bridges)
        lant.core = "MAGICHATCT"
        self.assertEqual(2, lant.number_of_lan_bridges)
        lant.core = "MAGICHATCCS"
        self.assertEqual(2, lant.number_of_lan_bridges)

    def test_monoisotopic_mass(self):
        "Test Lanthipeptide.monoisotopic_mass"
        lant = Lanthipeptide(23, 42, 17, 23, 'Class-I')
        lant.core = "MAGICHAT"
        analysis = ProteinAnalysis("MAGICHAT", monoisotopic=True)
        weight = analysis.molecular_weight()
        # Thr is assumed to be dehydrated
        weight -= 18
        self.assertAlmostEqual(weight, lant.monoisotopic_mass)
        self.assertAlmostEqual(weight, lant._monoisotopic_weight)

    def test_molecular_weight(self):
        "Test Lanthipeptide.molecular_weight"
        lant = Lanthipeptide(23, 42, 17, 23, 'Class-I')
        lant.core = "MAGICHAT"
        analysis = ProteinAnalysis("MAGICHAT", monoisotopic=False)
        weight = analysis.molecular_weight()
        # Thr is assumed to be dehydrated
        weight -= 18.02
        self.assertAlmostEqual(weight, lant.molecular_weight)
        self.assertAlmostEqual(weight, lant._weight)

    def test_alternative_weights(self):
        "Test Lanthipeptide.alt_weights"
        lant = Lanthipeptide(23, 42, 17, 23, 'Class-I')
        lant.core = "MAGICHATS"
        analysis = ProteinAnalysis("MAGICHATS", monoisotopic=False)
        weight = analysis.molecular_weight()
        # One Ser/Thr is assumed to be dehydrated, but not the other
        weight -= 18.02
        self.assertEqual([weight], lant.alternative_weights)


class TestSpecificAnalysis(unittest.TestCase):
    class FakeHit(object): # TODO: see antismash.common.test.helpers
        class FakeHsp(object): # TODO: see antismash.common.test.helpers
            def __init__(self, start, end, score):
                self.query_start = start
                self.query_end = end
                self.bitscore = score

        def __init__(self, start, end, score, desc):
            self.hsps = [self.FakeHsp(start, end, score)]
            self.description = desc

        def __iter__(self):
            return iter(self.hsps)

    def setUp(self):
        self.trace_tracker = TraceTracker()
        self.hmmpfam_return_vals = []
        mock('subprocessing.run_hmmpfam2', tracker=self.trace_tracker,
             returns=self.hmmpfam_return_vals)

    def tearDown(self):
        restore()

    def test_predict_cleavage_site(self):
        "Test lanthipeptides.predict_cleavage_site()"
        resvec = predict_cleavage_site('foo', 'bar')
        self.assertEqual(None, resvec)
        fake_hit = self.FakeHit(24, 42, 17, 'fake')
        self.hmmpfam_return_vals.append([fake_hit])

        res = predict_cleavage_site('foo', 'bar')
        self.assertEqual(23, res.start)
        self.assertEqual(42, res.end)
        self.assertEqual(17, res.score)
        self.assertEqual('fake', res.lantype)

    def test_result_vec_to_features(self):
        "Test lanthipeptides.result_vec_to_features()"
        orig_feature = DummyCDS(0, 165)
        orig_feature.locus_tag = 'FAKE0001'
        vec = Lanthipeptide(17, 23, 42, 23, 'Class-I')
        seq = "TAILTAILTAILTAILTAILTAILTAILTAILTAILCC"
        vec.core = seq
        vec.leader = "HEADHEADHEAD"
        new_feature = result_vec_to_feature(orig_feature, vec)
        new_features = new_feature.to_biopython()
        self.assertEqual(2, len(new_features))
        leader, core = new_features

        self.assertEqual(0, leader.location.start)
        self.assertEqual((23 * 3), leader.location.end)
        self.assertEqual(leader.location.strand, 1)
        self.assertEqual('CDS_motif', leader.type)
        self.assertEqual(orig_feature.locus_tag, leader.qualifiers['locus_tag'])
        self.assertEqual(['leader peptide', 'lanthipeptide',
                          'predicted leader seq: HEADHEADHEAD'], leader.qualifiers['note'])

        self.assertEqual(leader.location.end, core.location.start)
        self.assertEqual(165, core.location.end)
        self.assertEqual(1, core.location.strand)
        self.assertEqual('CDS_motif', core.type)
        expected = ['core peptide', 'lanthipeptide', 'monoisotopic mass: 3646.3',
                    'molecular weight: 3648.6',
                    'alternative weights: 3666.6; 3684.6; 3702.7; 3720.7; 3738.7; 3756.7; 3774.7',
                    'number of bridges: 2',
                    'predicted core seq: TAILTAILTAILTAILTAILTAILTAILTAILTAILCC',
                    'predicted class: Class I',
                    'score: 42.00',
                    'RODEO score: 23',
                   ]
        self.assertEqual(set(expected), set(core.qualifiers['note']))
        self.assertEqual(orig_feature.locus_tag, core.qualifiers['locus_tag'][0])
