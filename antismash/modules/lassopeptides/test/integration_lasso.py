# License: GNU Affero General Public License v3 or later
# A copy of GNU AGPL v3 should have been included in this software package in LICENSE.txt.

# for test files, silence irrelevant and noisy pylint warnings
# pylint: disable=no-self-use,protected-access,missing-docstring

import unittest

from helperlibs.bio import seqio

import antismash
from antismash.common import path, secmet
from antismash.common.test import helpers
from antismash.config import build_config, destroy_config
from antismash.modules import lassopeptides
from antismash.modules.lassopeptides import specific_analysis, check_prereqs

check_prereqs()  # ensure fimo detected


class IntegrationLasso(unittest.TestCase):
    def setUp(self):
        self.options = build_config(["--minimal", "--enable-lassopeptides"],
                                    isolated=True, modules=antismash.get_all_modules())

    def tearDown(self):
        destroy_config()

    def test_astexin1(self):
        "Test lassopeptides prediction for astexin-1"
        result = helpers.run_and_regenerate_results_for_module(path.get_full_path(__file__, "data", "astex1.gbk"),
                        lassopeptides, self.options, expected_record_count=1)
        assert list(result.motifs_by_locus) == ["ctg1_orf03094"]
        assert len(result.motifs_by_locus["ctg1_orf03094"]) == 1
        motif = result.motifs_by_locus["ctg1_orf03094"][0]
        # values considering tail cut
        self.assertAlmostEqual(motif.cut_mass, 2093.0, delta=1)
        self.assertAlmostEqual(motif.cut_weight, 2094.2, delta=1)
        assert motif.num_bridges == 0
        assert motif.leader == "MHTPIISETVQPKTAGLIVLGKASAETR"
        assert motif.core == "GLSQGVEPDIGQTYFEESR"
        assert motif.tail == "INQD"
        assert motif.macrolactam == 'GLSQGVEPD'
        assert motif.peptide_subclass == "Class II"
        self.assertAlmostEqual(motif.monoisotopic_mass, 2563.2, delta=1)
        self.assertAlmostEqual(motif.molecular_weight, 2564.7, delta=1)

    def test_burhizin(self):
        "Test lassopeptide prediction for burhizin"
        result = helpers.run_and_regenerate_results_for_module(path.get_full_path(__file__, "data", "burhizin.gbk"),
                        lassopeptides, self.options, expected_record_count=1)

        assert list(result.motifs_by_locus) == ["ctg1_orf02117"]
        assert len(result.motifs_by_locus["ctg1_orf02117"]) == 1
        motif = result.motifs_by_locus["ctg1_orf02117"][0]

        self.assertAlmostEqual(motif.cut_mass, 1846.9, delta=1)
        self.assertAlmostEqual(motif.cut_weight, 1847.9, delta=1)
        self.assertEqual(motif.num_bridges, 0)
        self.assertEqual(motif.leader, "MNKQQQESGLLLAEESLMELCASSETL")
        self.assertEqual(motif.core, "GGAGQYKEVEAGRWSDR")
        self.assertEqual(motif.peptide_subclass, 'Class II')
        self.assertEqual(motif.macrolactam, 'GGAGQYKE')
        self.assertEqual(motif.tail, 'IDSDDE')
        self.assertAlmostEqual(motif.monoisotopic_mass, 2521.1, delta=1)
        self.assertAlmostEqual(motif.molecular_weight, 2522.5, delta=1)

    def test_ssv2083(self):
        "Test lassopeptide prediction for SSV-2083"
        result = helpers.run_and_regenerate_results_for_module(path.get_full_path(__file__, "data", 'SSV2083.gbk'),
                        lassopeptides, self.options, expected_record_count=1)

        assert list(result.motifs_by_locus) == ["ctg1_orf11610"]
        assert len(result.motifs_by_locus["ctg1_orf11610"]) == 1
        motif = result.motifs_by_locus["ctg1_orf11610"][0]

        self.assertAlmostEqual(motif.cut_mass, 2086.8, delta=1)
        self.assertAlmostEqual(motif.cut_weight, 2088.4, delta=1)
        self.assertEqual(motif.num_bridges, 2)
        self.assertEqual(motif.leader, "VLISTTNGQGTPMTSTDELYEAPELIEIGDYAELTR")
        self.assertEqual(motif.core, "CVWGGDCTDFLGCGTAWICV")
        assert motif.macrolactam == "CVWGGDCTD"
        self.assertEqual(motif.tail, "")
        self.assertEqual(motif.peptide_subclass, 'Class I')
