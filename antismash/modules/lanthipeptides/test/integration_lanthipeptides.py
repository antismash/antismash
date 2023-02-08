# License: GNU Affero General Public License v3 or later
# A copy of GNU AGPL v3 should have been included in this software package in LICENSE.txt.

# for test files, silence irrelevant and noisy pylint warnings
# pylint: disable=use-implicit-booleaness-not-comparison,protected-access,missing-docstring

import os
import unittest

import antismash
from antismash.common import path
from antismash.common.test import helpers
from antismash.common.secmet import Record
from antismash.config import build_config, update_config, destroy_config
from antismash.modules import lanthipeptides
import antismash.modules.lanthipeptides.config as lanthi_config


class IntegrationLanthipeptides(unittest.TestCase):
    def setUp(self):
        self.options = build_config(["--minimal", "--enable-html", "--enable-lanthipeptides"],
                                    isolated=True, modules=antismash.get_all_modules())
        self.set_fimo_enabled(True)

    def tearDown(self):
        destroy_config()

    def set_fimo_enabled(self, val):
        update_config({"without_fimo": not val})
        lanthi_config.get_config().fimo_present = val

    def gather_all_motifs(self, result):
        motifs = []
        for locii in result.clusters.values():
            for locus in locii:
                motifs.extend(result.motifs_by_locus[locus])
        return motifs

    def run_lanthi(self, genbank, html_snippet, expected_motifs=1):
        def callback(output_dir):
            with open(os.path.join(output_dir, "index.html"), encoding="utf-8") as handle:
                content = handle.read()
                assert html_snippet in content

        rec = Record.from_genbank(genbank, taxon="bacteria")[0]
        assert not rec.get_cds_motifs()
        result = helpers.run_and_regenerate_results_for_module(genbank, lanthipeptides,
                                                               self.options, callback=callback)
        assert len(result.clusters) == 1
        motifs = self.gather_all_motifs(result)
        assert len(motifs) == expected_motifs

        result.add_to_record(rec)
        assert len(rec.get_cds_motifs()) == expected_motifs

        return rec, result

    def test_nisin(self):
        "Test lanthipeptide prediction for nisin A"
        expected_html_snippet = "Leader:"
        genbank = helpers.get_path_to_nisin_genbank()
        rec, _ = self.run_lanthi(genbank, expected_html_snippet)

        prepeptide = rec.get_cds_motifs()[0]
        # real monoisotopic mass is 3351.51, but we overpredict a Dha
        self.assertAlmostEqual(3333.6, prepeptide.monoisotopic_mass, delta=0.05)
        # real mw is 3354.5, see above
        self.assertAlmostEqual(3336.0, prepeptide.molecular_weight, delta=0.05)
        for expected, calculated in zip([3354.0, 3372.1, 3390.1, 3408.1],
                                        prepeptide.alternative_weights):
            self.assertAlmostEqual(expected, calculated, delta=0.05)
        assert prepeptide.detailed_information.lan_bridges == 5
        self.assertEqual("MSTKDFNLDLVSVSKKDSGASPR", prepeptide.leader)
        self.assertEqual("ITSISLCTPGCKTGALMGCNMKTATCHCSIHVSK", prepeptide.core)
        self.assertEqual('Class I', prepeptide.peptide_subclass)

    def test_epidermin(self):
        "Test lanthipeptide prediction for epidermin"
        filename = path.get_full_path(__file__, 'data', 'epidermin.gbk')
        expected_html_snippet = "Lanthipeptide(s) for epiB"
        rec, _ = self.run_lanthi(filename, expected_html_snippet)

        prepeptide = rec.get_cds_motifs()[0]
        self.assertAlmostEqual(2164, prepeptide.monoisotopic_mass, delta=0.5)
        self.assertAlmostEqual(2165.6, prepeptide.molecular_weight, delta=0.5)
        self.assertEqual(3, prepeptide.detailed_information.lan_bridges)
        self.assertEqual("MEAVKEKNDLFNLDVKVNAKESNDSGAEPR", prepeptide.leader)
        self.assertEqual("IASKFICTPGCAKTGSFNSYCC", prepeptide.core)
        self.assertEqual('Class I', prepeptide.peptide_subclass)
        self.assertEqual(['AviCys'], prepeptide.detailed_information.get_modifications())

    def test_microbisporicin(self):
        "Test lanthipeptide prediction for microbisporicin"
        filename = path.get_full_path(__file__, 'data', 'microbisporicin.gbk')
        expected_snippet = "Lanthipeptide(s) for mibB"
        rec, _ = self.run_lanthi(filename, expected_snippet)

        prepeptide = rec.get_cds_motifs()[0]
        # NOTE: this is not the correct weight for microbisporicin
        # there are some additional modifications we do not predict yet
        self.assertAlmostEqual(2212.9, prepeptide.monoisotopic_mass, delta=0.5)
        self.assertAlmostEqual(2214.5, prepeptide.molecular_weight, delta=0.5)
        self.assertEqual(4, prepeptide.detailed_information.lan_bridges)
        self.assertEqual("MPADILETRTSETEDLLDLDLSIGVEEITAGPA", prepeptide.leader)
        self.assertEqual("VTSWSLCTPGCTSPGGGSNCSFCC", prepeptide.core)
        self.assertEqual('Class I', prepeptide.peptide_subclass)
        self.assertEqual(['AviCys', 'Cl', 'OH'], prepeptide.detailed_information.get_modifications())

    def test_epicidin(self):
        "Test lanthipeptide prediction for epicidin 280"
        filename = path.get_full_path(__file__, 'data', 'epicidin_280.gbk')
        expected_html_snippet = "Core:"
        _, result = self.run_lanthi(filename, expected_html_snippet)
        prepeptide = self.gather_all_motifs(result)[0]
        self.assertAlmostEqual(3115.7, prepeptide.monoisotopic_mass, delta=0.5)
        self.assertAlmostEqual(3117.7, prepeptide.molecular_weight, delta=0.5)
        for expected, calculated in zip([3135.7, 3153.7, 3171.7],
                                        prepeptide.alternative_weights):
            self.assertAlmostEqual(expected, calculated, delta=0.05)
        self.assertEqual(3, prepeptide.detailed_information.lan_bridges)
        self.assertEqual("MENKKDLFDLEIKKDNMENNNELEAQ", prepeptide.leader)
        self.assertEqual("SLGPAIKATRQVCPKATRFVTVSCKKSDCQ", prepeptide.core)
        self.assertEqual('Class I', prepeptide.peptide_subclass)
        self.assertEqual(['Lac'], prepeptide.detailed_information.get_modifications())

    def test_labyrinthopeptin(self):
        "Test lanthipeptide prediction for labyrinthopeptin"
        filename = path.get_full_path(__file__, 'data', 'labyrinthopeptin.gbk')
        expected_snippet = "Core:"
        rec, _ = self.run_lanthi(filename, expected_snippet, expected_motifs=2)
        motifs = sorted(rec.get_cds_motifs(), key=lambda x: x.locus_tag)

        assert motifs[0].locus_tag == "labA1_lanthipeptide"
        assert motifs[0].detailed_information.lan_bridges == 4
        self.assertAlmostEqual(motifs[0].molecular_weight, 2059.3, delta=0.5)

        assert motifs[1].locus_tag == "labA2_lanthipeptide"
        assert motifs[1].detailed_information.lan_bridges == 4
        self.assertAlmostEqual(motifs[1].molecular_weight, 1908.1, delta=0.5)

    def test_sco_cluster3(self):
        "Test lanthipeptide prediction for SCO cluster #3"
        filename = path.get_full_path(__file__, 'data', 'sco_cluster3.gbk')
        expected_snippet = "Core:"
        rec, _ = self.run_lanthi(filename, expected_snippet)

        prepeptide = rec.get_cds_motifs()[0]
        assert prepeptide.peptide_subclass == 'Class I'
        assert prepeptide.detailed_information.lan_bridges == 3
        assert prepeptide.detailed_information.rodeo_score == 18

    def test_lactocin_s(self):
        """Test lanthipeptide prediction for lactocin S"""
        filename = path.get_full_path(__file__, 'data', 'lactocin_s.gbk')
        expected_snippet = "Core:"
        _, result = self.run_lanthi(filename, expected_snippet)

        assert len(result.clusters) == 1
        assert result.clusters[1] == set(["lasM"])
        motifs = result.motifs_by_locus["lasM"]
        assert len(motifs) == 1

        assert motifs[0].peptide_subclass == "Class II"
        assert motifs[0].locus_tag == "lasA_lanthipeptide"

    def test_multiple_biosynthetic_enzymes(self):
        # TODO: find/create an input with both class II and class III lanthipeptides
        # this was the case in CP013129.1, in a nrps-lanthipeptide hybrid, but
        # the hybrid was only created due to a bug in cluster formation
        pass


class IntegrationLanthipeptidesWithoutFimo(IntegrationLanthipeptides):
    def setUp(self):
        self.options = build_config(["--minimal", "--enable-html", "--enable-lanthipeptides"],
                                    isolated=True, modules=antismash.get_all_modules())
        self.set_fimo_enabled(False)
        assert lanthi_config.get_config().fimo_present is False
