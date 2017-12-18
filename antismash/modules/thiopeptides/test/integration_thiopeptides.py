# License: GNU Affero General Public License v3 or later
# A copy of GNU AGPL v3 should have been included in this software package in LICENSE.txt.

# for test files, silence irrelevant and noisy pylint warnings
# pylint: disable=no-self-use,protected-access,missing-docstring

import os
import unittest

from helperlibs.bio import seqio

import antismash
from antismash.common import path, secmet
from antismash.common.test import helpers
from antismash.config import update_config, destroy_config, build_config
from antismash.modules import thiopeptides


class TestIntegration(unittest.TestCase):
    def setUp(self):
        update_config({"cpus": 1})

    def tearDown(self):
        destroy_config()

    def test_nosiheptide(self):
        "Test thiopeptide prediction for nosiheptide - nosM"
        rec = seqio.read(path.get_full_path(__file__, 'data', 'nosi_before_analysis.gbk'))
        rec = secmet.Record.from_biopython(rec)
        rec.get_cluster(1).trim_overlapping()
        assert rec.get_feature_count() == 56
        assert not rec.get_cds_motifs()
        result = thiopeptides.specific_analysis(rec)
        assert rec.get_feature_count() == 56

        assert len(result.motifs) == 1

        result.add_to_record(rec)
        for i in rec.get_cds_motifs():
            print(i, i.leader, i.score, i.rodeo_score)
        assert len(rec.get_cds_motifs()) == 1, rec.get_cds_motifs()
        assert rec.get_feature_count() == 57

        # check the motif in an existing CDS
        prepeptide = rec.get_cds_motifs()[0]
        assert prepeptide is result.motifs[0]

        self.assertAlmostEqual(1315.3, prepeptide.monoisotopic_mass, places=1)
        self.assertAlmostEqual(1316.5, prepeptide.molecular_weight, places=1)
        assert prepeptide.leader == "MDAAHLSDLDIDALEISEFLDESRLEDSEVVAKVMSA"
        assert prepeptide.core == "SCTTCECCCSCSS"
        assert prepeptide.macrocycle == "26-member"
        assert prepeptide.peptide_subclass == "Type I"
        self.assertAlmostEqual(1222.4, prepeptide.mature_weights[0], places=1)
        self.assertAlmostEqual(1221.2, prepeptide.mature_weights[1], places=1)
        for calc, expected in zip(prepeptide.mature_weights[2:],
                                  [1240.4, 1258.4, 1276.5, 1294.5, 1312.5, 1330.5]):
            self.assertAlmostEqual(calc, expected, places=1)
        assert prepeptide.core_features == 'Central ring: pyridine tetrasubstituted (hydroxyl group present); second macrocycle'
        assert prepeptide.tail_reaction == 'dealkylation of C-Terminal residue; amidation'

    def test_lactazole(self):
        "Test thiopeptide prediction for lactazole - lazA"
        rec = seqio.read(path.get_full_path(__file__, 'data', 'lac_before_analysis.gbk'))
        rec = secmet.Record.from_biopython(rec)
        assert rec.get_feature_count() == 21
        assert not rec.get_cds_motifs()
        results = thiopeptides.specific_analysis(rec)
        assert len(results.motifs) == 1
        # ensure record not adjusted yet
        assert rec.get_feature_count() == 21
        assert not rec.get_cds_motifs()

        # add and check new motif added
        results.add_to_record(rec)
        assert rec.get_feature_count() == 22
        assert len(rec.get_cds_motifs()) == 1
        prepeptide = rec.get_cds_motifs()[0]
        assert prepeptide is results.motifs[0]
        self.assertAlmostEqual(1362.5, prepeptide.monoisotopic_mass, places=1)
        self.assertAlmostEqual(1363.5, prepeptide.molecular_weight, places=1)
        assert prepeptide.leader == "MSDITASRVESLDLQDLDLSELTVTSLRDTVALPENGA"
        assert prepeptide.core == "SWGSCSCQASSSCA"
        assert not prepeptide.macrocycle
        assert prepeptide.peptide_subclass == "Type III"
        assert prepeptide.core_features == 'Central ring: pyridine trisubstituted'
        assert prepeptide.tail == 'QPQDM'
        for calc, expected in zip(prepeptide.alternative_weights,
                                  [1381.5, 1399.5, 1417.5, 1435.5, 1453.6, 1471.6]):
            self.assertAlmostEqual(calc, expected, places=1)
        assert len(prepeptide.to_biopython()) == 3  # leader, core, tail

    def test_thiostrepton(self):
        "Test thiopeptide prediction for thiostrepton"
        rec = seqio.read(path.get_full_path(__file__, 'data', 'thiostrepton_before_analysis.gbk'))
        rec = secmet.Record.from_biopython(rec)
        assert rec.get_feature_count() == 27
        # two existing motifs
        assert len(rec.get_cds_motifs()) == 2

        results = thiopeptides.specific_analysis(rec)
        assert len(results.motifs) == 1
        # ensure record not adjusted yet
        self.assertEqual(27, rec.get_feature_count())
        results.add_to_record(rec)
        # the new motif is added
        assert rec.get_feature_count() == 28
        assert len(rec.get_cds_motifs()) == 3
        prepeptides = rec.get_cds_motifs()
        prepeptide = None
        for feature in prepeptides:
            if isinstance(feature, secmet.feature.Prepeptide):
                prepeptide = feature
                break
        # and the motif that was added is exactly the one in results
        assert prepeptide is results.motifs[0]
        self.check_thiostrepton_values(prepeptide)

    def check_thiostrepton_values(self, prepeptide):
        # check the values are as expected
        self.assertAlmostEqual(1639.6, prepeptide.monoisotopic_mass, places=1)
        self.assertAlmostEqual(1640.9, prepeptide.molecular_weight, places=1)
        assert prepeptide.leader == "MSNAALEIGVEGLTGLDVDTLEISDYMDETLLDGEDLTVTM"
        assert prepeptide.core == "IASASCTTCICTCSCSS"
        assert prepeptide.macrocycle == "26-member"
        assert prepeptide.peptide_subclass == "Type II"
        assert prepeptide.core_features == 'Central ring: piperidine; second macrocycle containing a quinaldic acid moiety'
        self.assertAlmostEqual(1646.8, prepeptide.mature_weights[0], places=1)
        self.assertAlmostEqual(1645.5, prepeptide.mature_weights[1], places=1)
        for calc, expect in zip(prepeptide.mature_weights[2:],
                                [1664.8, 1682.8, 1700.8, 1718.8, 1736.9, 1754.9, 1772.9, 1790.9]):
            self.assertAlmostEqual(calc, expect, places=1)

    def test_thiostrepton_full(self):
        def callback(outdir):
            # make sure the html_output section was tested
            with open(os.path.join(outdir, "index.html")) as handle:
                content = handle.read()
                assert "ACN52291.1 leader / core peptide, putative Type II" in content
        config = build_config(["--minimal", "--enable-thiopeptides"],
                              isolated=True, modules=antismash.get_all_modules())

        record_path = path.get_full_path(__file__, 'data', 'thiostrepton_before_analysis.gbk')
        regenned_results = helpers.run_and_regenerate_results_for_module(record_path, thiopeptides, config,
                                     expected_record_count=1, callback=callback)
        assert regenned_results
        assert len(regenned_results.motifs) == 1
        self.check_thiostrepton_values(regenned_results.motifs[0])

    def test_CP009369(self):  # pylint: disable=invalid-name
        " tests the special case HMM files for rodeo "
        config = build_config(["--minimal", "--enable-thiopeptides"],
                              isolated=True, modules=antismash.get_all_modules())
        record_path = path.get_full_path(__file__, 'data', 'CP009369.1.gbk')
        results = helpers.run_and_regenerate_results_for_module(record_path, thiopeptides, config,
                                     expected_record_count=1)
        assert results
        assert len(results.motifs) == 1
        prepeptide = results.motifs[0]
        self.assertAlmostEqual(1934.6, prepeptide.monoisotopic_mass, places=1)
        self.assertAlmostEqual(1936.0, prepeptide.molecular_weight, places=1)
        assert prepeptide.leader == "MVKSIIKARESGRFYETKYLKGGEEMKEQKELKNEEFELDVEFLDLDEVSAIPETTA"
        assert prepeptide.core == "SSGTSSCSASSTCGSSSCCGSC"
        assert not prepeptide.macrocycle
        assert prepeptide.peptide_subclass == "Type III"
        assert prepeptide.core_features == 'Central ring: pyridine trisubstituted'
        assert prepeptide.tail == ''
        for calc, expected in zip(prepeptide.alternative_weights,
                                  [1954.0, 1972.1, 1990.1, 2008.1, 2026.1, 2044.1,
                                   2062.2, 2080.2, 2098.2, 2116.2, 2134.2, 2152.3, 2170.3]):
            self.assertAlmostEqual(calc, expected, places=1)
        assert len(prepeptide.to_biopython()) == 2  # no tail
