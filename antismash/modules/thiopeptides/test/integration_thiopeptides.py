# License: GNU Affero General Public License v3 or later
# A copy of GNU AGPL v3 should have been included in this software package in LICENSE.txt.

# for test files, silence irrelevant and noisy pylint warnings
# pylint: disable=use-implicit-booleaness-not-comparison,protected-access,missing-docstring

import os
import unittest

import antismash
from antismash.common import path, secmet
from antismash.common.test import helpers
from antismash.config import destroy_config, build_config
from antismash.modules import thiopeptides


class TestIntegration(unittest.TestCase):
    def setUp(self):
        self.config = build_config(["--minimal", "--enable-html", "--enable-thiopeptides", "--cpus", "1"],
                                   isolated=True, modules=antismash.get_all_modules())

    def tearDown(self):
        destroy_config()

    def test_nosiheptide(self):
        "Test thiopeptide prediction for nosiheptide - nosM"
        genbank = path.get_full_path(__file__, 'data', 'nosi_before_analysis.gbk')
        result = helpers.run_and_regenerate_results_for_module(genbank, thiopeptides, self.config)
        rec = secmet.Record.from_genbank(genbank, "bacteria")[0]
        existing_feature_count = rec.get_feature_count()
        assert result.motifs[0].peptide_subclass == "Type I"
        assert not rec.get_cds_motifs()
        assert rec.get_feature_count() == existing_feature_count

        assert len(result.motifs) == 1

        result.add_to_record(rec)
        for i in rec.get_cds_motifs():
            print(i, i.leader, i.score, i.detailed_information.rodeo_score)
        assert len(rec.get_cds_motifs()) == 1, rec.get_cds_motifs()
        assert rec.get_feature_count() == existing_feature_count + 1

        # check the motif in an existing CDS
        prepeptide = rec.get_cds_motifs()[0]
        assert prepeptide is result.motifs[0]

        self.assertAlmostEqual(1315.3, prepeptide.monoisotopic_mass, places=1)
        self.assertAlmostEqual(1316.5, prepeptide.molecular_weight, places=1)
        assert prepeptide.leader == "MDAAHLSDLDIDALEISEFLDESRLEDSEVVAKVMSA"
        assert prepeptide.core == "SCTTCECCCSCSS"
        assert prepeptide.detailed_information.macrocycle == "26-member"
        assert prepeptide.peptide_subclass == "Type I"
        for calc, expected in zip(prepeptide.detailed_information.mature_weights,
                                  [1222.4, 1221.2, 1240.4, 1258.4, 1276.5, 1294.5, 1312.5, 1330.5]):
            self.assertAlmostEqual(calc, expected, places=1)
        expected_core_features = ("Central ring: pyridine tetrasubstituted (hydroxyl group present);"
                                  " second macrocycle")
        assert prepeptide.detailed_information.core_features == expected_core_features
        assert prepeptide.detailed_information.amidation

    def test_lactazole(self):
        "Test thiopeptide prediction for lactazole - lazA"
        genbank = path.get_full_path(__file__, 'data', 'lac_before_analysis.gbk')

        results = helpers.run_and_regenerate_results_for_module(genbank, thiopeptides, self.config)
        assert len(results.motifs) == 1

        rec = secmet.Record.from_genbank(genbank, "bacteria")[0]
        before_count = rec.get_feature_count()
        assert not rec.get_cds_motifs()

        # add and check new motif added
        results.add_to_record(rec)
        assert rec.get_feature_count() == before_count + 1
        assert len(rec.get_cds_motifs()) == 1
        prepeptide = rec.get_cds_motifs()[0]
        assert prepeptide is results.motifs[0]

        # check prepeptide values
        self.assertAlmostEqual(1362.5, prepeptide.monoisotopic_mass, places=1)
        self.assertAlmostEqual(1363.5, prepeptide.molecular_weight, places=1)
        assert prepeptide.leader == "MSDITASRVESLDLQDLDLSELTVTSLRDTVALPENGA"
        assert prepeptide.core == "SWGSCSCQASSSCA"
        assert not prepeptide.detailed_information.macrocycle
        assert prepeptide.peptide_subclass == "Type III"
        assert prepeptide.detailed_information.core_features == 'Central ring: pyridine trisubstituted'
        assert prepeptide.tail == 'QPQDM'
        for calc, expected in zip(prepeptide.alternative_weights,
                                  [1381.5, 1399.5, 1417.5, 1435.5, 1453.6, 1471.6]):
            self.assertAlmostEqual(calc, expected, places=1)
        assert len(prepeptide.to_biopython()) == 3  # leader, core, tail

    def test_thiostrepton(self):
        def callback(outdir):
            # make sure the html_output section was tested
            with open(os.path.join(outdir, "index.html"), encoding="utf-8") as handle:
                content = handle.read()
                assert "Leader:" in content

        genbank = path.get_full_path(__file__, 'data', 'thiostrepton_before_analysis.gbk')
        result = helpers.run_and_regenerate_results_for_module(genbank, thiopeptides,
                                                               self.config, callback=callback)

        result = helpers.run_and_regenerate_results_for_module(genbank, thiopeptides, self.config)

        assert len(result.motifs) == 1
        prepeptide = result.motifs[0]

        # check prepeptide values
        self.assertAlmostEqual(1639.6, prepeptide.monoisotopic_mass, places=1)
        self.assertAlmostEqual(1640.9, prepeptide.molecular_weight, places=1)
        assert prepeptide.leader == "MSNAALEIGVEGLTGLDVDTLEISDYMDETLLDGEDLTVTM"
        assert prepeptide.core == "IASASCTTCICTCSCSS"
        assert prepeptide.detailed_information.macrocycle == "26-member"
        assert prepeptide.peptide_subclass == "Type II"
        expected_features = ("Central ring: piperidine;"
                             " second macrocycle containing a quinaldic acid moiety")
        assert prepeptide.detailed_information.core_features == expected_features
        for calc, expect in zip(prepeptide.detailed_information.mature_weights,
                                [1646.8, 1645.5, 1664.8, 1682.8, 1700.8, 1718.8,
                                 1736.9, 1754.9, 1772.9, 1790.9]):
            self.assertAlmostEqual(calc, expect, places=1)

    def test_CP009369(self):  # pylint: disable=invalid-name
        " tests the special case HMM files for rodeo "
        record_path = path.get_full_path(__file__, 'data', 'CP009369.1.gbk')
        results = helpers.run_and_regenerate_results_for_module(record_path, thiopeptides, self.config)
        assert results
        assert len(results.motifs) == 2
        prepeptide = results.motifs[0]
        other = results.motifs[1]
        self.assertAlmostEqual(1934.6, prepeptide.monoisotopic_mass, places=1)
        self.assertAlmostEqual(1936.0, prepeptide.molecular_weight, places=1)
        assert prepeptide.leader == "MVKSIIKARESGRFYETKYLKGGEEMKEQKELKNEEFELDVEFLDLDEVSAIPETTA"
        assert other.leader == "MDVEFLDLDEVSAIPETTA"
        assert prepeptide.core == "SSGTSSCSASSTCGSSSCCGSC"
        assert other.core == prepeptide.core
        assert not prepeptide.detailed_information.macrocycle
        assert prepeptide.peptide_subclass == "Type III"
        assert prepeptide.detailed_information.core_features == 'Central ring: pyridine trisubstituted'
        assert prepeptide.tail == ''
        for calc, expected in zip(prepeptide.alternative_weights,
                                  [1954.0, 1972.1, 1990.1, 2008.1, 2026.1, 2044.1,
                                   2062.2, 2080.2, 2098.2, 2116.2, 2134.2, 2152.3, 2170.3]):
            self.assertAlmostEqual(calc, expected, places=1)
        assert len(prepeptide.to_biopython()) == 2  # no tail

    def test_NZ_CP015439(self):  # pylint: disable=invalid-name
        """ Tests that small ORFs are found and saved in results """
        record_path = path.get_full_path(__file__, 'data', 'NZ_CP015439_section.gbk')
        results = helpers.run_and_regenerate_results_for_module(record_path, thiopeptides, self.config)
        assert results

        # check that the extra orf was found and stored correctly
        assert len(results._cds_features) == 1
        additions = list(results._cds_features.values())[0]
        assert len(additions) == 1
        assert isinstance(additions[0], secmet.features.CDSFeature)

        # also test the analysis results itself
        assert len(results.motifs) == 1
        prepeptide = results.motifs[0]
        self.assertAlmostEqual(1408.6, prepeptide.monoisotopic_mass, places=1)
        self.assertAlmostEqual(1409.5, prepeptide.molecular_weight, places=1)
        assert prepeptide.leader == "MPDITQYTAAGTSTLSTESSVLSASCP"
        assert prepeptide.core == "TSTASTYTSMSSVS"
