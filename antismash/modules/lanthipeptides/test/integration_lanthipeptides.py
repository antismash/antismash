# License: GNU Affero General Public License v3 or later
# A copy of GNU AGPL v3 should have been included in this software package in LICENSE.txt.

import os
import unittest

from helperlibs.bio import seqio
from helperlibs.wrappers.io import TemporaryDirectory

import antismash
from antismash.common import path
from antismash.common.test import helpers
from antismash.common.secmet import Record
from antismash.config.args import Config, build_parser
from antismash.modules import lanthipeptides
from antismash.modules.lanthipeptides import specific_analysis, LanthiResults

class IntegrationLantipeptides(unittest.TestCase):
    def setUp(self):
        options = helpers.get_simple_options(lanthipeptides, [])
        options.without_fimo = False
        Config(options)

    def tearDown(self):
        Config({})

    def test_nisin(self):
        "Test lanthipeptide prediction for nisin A"
        rec = Record.from_biopython(seqio.read(helpers.get_path_to_nisin_with_detection()))
        assert not rec.get_cds_motifs()
        result = specific_analysis(rec)
        assert len(result.motifs) == 1
        assert len(result.clusters_with_motifs) == 1
        assert len(rec.get_cds_motifs()) == 1
        prepeptide = result.motifs[0]
        # real monoisotopic mass is 3351.51, but we overpredict a Dha
        self.assertAlmostEqual(3333.6, prepeptide.monoisotopic_mass, delta=0.05)
        # real mw is 3354.5, see above
        self.assertAlmostEqual(3336.0, prepeptide.molecular_weight, delta=0.05)
        for expected, calculated in zip([3354.0, 3372.1, 3390.1, 3408.1],
                                        prepeptide.alternative_weights):
            self.assertAlmostEqual(expected, calculated, delta=0.05)
        assert prepeptide.lan_bridges == 5
        self.assertEqual("MSTKDFNLDLVSVSKKDSGASPR", prepeptide.leader_seq)
        self.assertEqual("ITSISLCTPGCKTGALMGCNMKTATCHCSIHVSK", prepeptide.core_seq)
        self.assertEqual('Class I', prepeptide.peptide_class)

        initial_json = result.to_json()
        regenerated = LanthiResults.from_json(initial_json, rec)
        assert len(result.motifs) == len(regenerated.motifs)
        assert len(result.clusters_with_motifs) == len(regenerated.clusters_with_motifs)
        assert result.motifs[0].location == regenerated.motifs[0].location
        assert initial_json == regenerated.to_json()

    def test_nisin_complete(self):
        with TemporaryDirectory() as output_dir:
            args = ["run_antismash.py", "--minimal", "--enable-lanthipeptides", "--output-dir", output_dir]
            parser = build_parser(modules=antismash.get_all_modules())
            options = parser.parse_args(args)
            options.all_enabled_modules = [module for module in antismash.get_all_modules() if module.is_enabled(options)] #TODO: shift elsewhere
            antismash.run_antismash(helpers.get_path_to_nisin_genbank(), Config(options))

            # make sure the html_output section was tested
            with open(os.path.join(output_dir, "index.html")) as handle:
                content = handle.read()
                assert "nisA leader / core peptide, putative Class I" in content


    def test_epidermin(self):
        "Test lanthipeptide prediction for epidermin"
        rec = Record.from_biopython(seqio.read(path.get_full_path(__file__, 'data/epidermin.gbk')))
        assert not rec.get_cds_motifs()
        result = specific_analysis(rec)
        assert len(rec.get_cds_motifs()) == 1
        assert len(result.motifs) == 1
        prepeptide = result.motifs[0]
        self.assertAlmostEqual(2164, prepeptide.monoisotopic_mass, delta=0.5)
        self.assertAlmostEqual(2165.6, prepeptide.molecular_weight, delta=0.5)
        self.assertEqual(3, prepeptide.lan_bridges)
        self.assertEqual("MEAVKEKNDLFNLDVKVNAKESNDSGAEPR", prepeptide.leader_seq)
        self.assertEqual("IASKFICTPGCAKTGSFNSYCC", prepeptide.core_seq)
        self.assertEqual('Class I', prepeptide.peptide_class)
        self.assertEqual(['AviCys'], prepeptide.get_modifications())

    def test_microbisporicin(self):
        "Test lanthipeptide prediction for microbisporicin"
        rec = Record.from_biopython(seqio.read(path.get_full_path(__file__, 'data/microbisporicin.gbk')))
        assert not rec.get_cds_motifs()
        result = specific_analysis(rec)
        assert len(result.motifs) == 1
        assert len(rec.get_cds_motifs()) == 1

        prepeptide = result.motifs[0]
        # NOTE: this is not the correct weight for microbisporicin
        # there are some additional modifications we do not predict yet
        self.assertAlmostEqual(2212.9, prepeptide.monoisotopic_mass, delta=0.5)
        self.assertAlmostEqual(2214.5, prepeptide.molecular_weight, delta=0.5)
        self.assertEqual(4, prepeptide.lan_bridges)
        self.assertEqual("MPADILETRTSETEDLLDLDLSIGVEEITAGPA", prepeptide.leader_seq)
        self.assertEqual("VTSWSLCTPGCTSPGGGSNCSFCC", prepeptide.core_seq)
        self.assertEqual('Class I', prepeptide.peptide_class)
        self.assertEqual(['AviCys', 'Cl', 'OH'], prepeptide.get_modifications())

    def test_epicidin(self):
        "Test lanthipeptide prediction for epicidin 280"
        rec = Record.from_biopython(seqio.read(path.get_full_path(__file__, 'data/epicidin_280.gbk')))
        assert not rec.get_cds_motifs()
        result = specific_analysis(rec)
        assert len(result.motifs) == 1
        assert len(rec.get_cds_motifs()) == 1

        prepeptide = result.motifs[0]
        self.assertAlmostEqual(3115.7, prepeptide.monoisotopic_mass, delta=0.5)
        self.assertAlmostEqual(3117.7, prepeptide.molecular_weight, delta=0.5)
        for expected, calculated in zip([3135.7, 3153.7, 3171.7],
                                        prepeptide.alternative_weights):
            self.assertAlmostEqual(expected, calculated, delta=0.05)
        self.assertEqual(3, prepeptide.lan_bridges)
        self.assertEqual("MENKKDLFDLEIKKDNMENNNELEAQ", prepeptide.leader_seq)
        self.assertEqual("SLGPAIKATRQVCPKATRFVTVSCKKSDCQ", prepeptide.core_seq)
        self.assertEqual('Class I', prepeptide.peptide_class)
        self.assertEqual(['Lac'], prepeptide.get_modifications())

    def test_labyrinthopeptin(self):
        "Test lanthipeptide prediction for labyrinthopeptin"
        rec = Record.from_biopython(seqio.read(path.get_full_path(__file__, 'data/labyrinthopeptin.gbk')))
        assert not rec.get_cds_motifs()
        result = specific_analysis(rec)
        assert len(result.motifs) == 2
        assert len(rec.get_cds_motifs()) == 2

    def test_sco_cluster3(self):
        "Test lanthipeptide prediction for SCO cluster #3"
        rec = Record.from_biopython(seqio.read(path.get_full_path(__file__, 'data/sco_cluster3.gbk')))
        assert not rec.get_cds_motifs()
        result = specific_analysis(rec)
        assert len(result.motifs) == 1
        assert len(rec.get_cds_motifs()) == 1
        self.assertEqual('Class I', result.motifs[0].peptide_class)

    def test_lactocin_s(self):
        """Test lanthipeptide prediction for lactocin S"""
        rec = Record.from_biopython(seqio.read(path.get_full_path(__file__, 'data/lactocin_s.gbk')))
        assert not rec.get_cds_motifs()
        result = specific_analysis(rec)
        assert len(result.motifs) == 1
        assert len(rec.get_cds_motifs()) == 1
        self.assertEqual('Class II', result.motifs[0].peptide_class)
