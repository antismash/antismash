# License: GNU Affero General Public License v3 or later
# A copy of GNU AGPL v3 should have been included in this software package in LICENSE.txt.

# for test files, silence irrelevant and noisy pylint warnings
# pylint: disable=use-implicit-booleaness-not-comparison,protected-access,missing-docstring

import unittest

import antismash
from antismash.common import fasta, path, secmet
from antismash.common.test import helpers
from antismash.config import build_config, destroy_config, get_config, update_config
from antismash.detection.genefunctions.tools.halogenases import (
    classify as halo_analysis,
    prepare_data,
)
from antismash.detection.genefunctions.tools.halogenases.data_structures import (
    Conventionality,
)

TRANSLATIONS = fasta.read_fasta(path.get_full_path(__file__, "data", "translations.fasta"))
# convert to just gene name, instead of the full info in each FASTA header
TRANSLATIONS = {key.rsplit("|", 1)[-1]: value for key, value in TRANSLATIONS.items()}


class TestHalongenases(unittest.TestCase):
    def setUp(self):
        options = build_config(self.get_args(), isolated=True,
                               modules=antismash.get_all_modules())
        self.old_config = get_config().__dict__
        self.options = update_config(options)
        prepare_data()

    def tearDown(self):
        destroy_config()
        update_config(self.old_config)

    def get_args(self):
        return ["--minimal", "--enable-genefunctions"]

    def generate_record(self, gene_name):
        features = [helpers.DummyCDS(locus_tag=gene_name, translation=TRANSLATIONS[gene_name])]
        record = helpers.DummyRecord(features=features)
        record.get_cds_features_within_regions = lambda: features
        return record

    def get_single_result(self, gene_name):
        record = self.generate_record(gene_name)
        results = halo_analysis(record.get_cds_features_within_regions(), self.options)

        # single best hits need to be present
        best_hits = results.best_hits
        assert len(best_hits) == 1
        hit = best_hits[gene_name]
        assert hit
        assert hit.subfunctions == ["Halogenation"]

        # and when present, function mapping must also exist and be of the same supertype
        assert results.function_mapping[gene_name] == secmet.qualifiers.GeneFunction.ADDITIONAL
        return hit

    def test_multiple(self):
        bha_name = "bhaA"
        chl_name = "ChlB4"
        record = secmet.Record.from_genbank(helpers.get_path_to_balhymicin_genbank())[0]
        # the real halogenase had best be present
        assert record.get_cds_by_name(bha_name)
        # and add a dummy second one, just to be sure multiple in a record will function correctly
        dummy_loc = secmet.locations.FeatureLocation(10, 10 + len(TRANSLATIONS["ChlB4"]) * 3, 1)
        dummy_halogenase = secmet.CDSFeature(dummy_loc, locus_tag=chl_name, translation=TRANSLATIONS[chl_name])
        record.add_cds_feature(dummy_halogenase)
        record.add_subregion(secmet.SubRegion(secmet.locations.FeatureLocation(0, len(record)), tool="dummy"))
        record.create_regions()
        assert dummy_halogenase in record.get_cds_features_within_regions()
        results = halo_analysis(record.get_cds_features_within_regions(), self.options)

        assert len(results.best_hits) == 2
        bha = results.best_hits[bha_name]
        chl_best = results.best_hits[chl_name]

        assert bha.query_id == bha_name
        assert len(bha.potential_matches) == 1
        match = bha.potential_matches[0]
        assert match.profile == "tyrosine-like_hpg_FDH"
        assert match.confidence == 1.0
        assert match.substrate == "Tyr"

        # ChlB4 is promiscuous, with two perfect matches against motifs from different groups of profiles
        # so the "best" hit should be the base conventional hit, with no substrates
        assert chl_best.query_id == chl_name
        assert chl_best.is_conventional()
        assert len(chl_best.potential_matches) == 0

        # however the substrate profile matches should also be included in all the hits
        chl_all = results.all_hits[chl_name]
        print([p.profile.name for p in chl_all])
        assert len(chl_all) == 3
        assert chl_all[0] is chl_best
        assert chl_all[1].confidence == chl_all[2].confidence
        assert {hit.profile.name for hit in chl_all[1:]} == {"orsellinic", "pyrrole"}

    def test_unconventional(self):
        result = self.get_single_result("VatD")
        assert not result.potential_matches
        assert not result.conventionality_residues
        assert result.conventionality == Conventionality.UNCONVENTIONAL
        assert result.confidence == 0.

    def test_generic_conventional(self):
        result = self.get_single_result("CtoA")
        assert not result.potential_matches
        assert result.conventionality_residues == {"W.W.I.": "WIWVIR"}
        assert result.conventionality == Conventionality.CONVENTIONAL
        assert result.confidence == 1.

    def test_tryptophan_7(self):
        result = self.get_single_result("ktzQ")

        assert len(result.potential_matches) == 1
        assert result.conventionality == Conventionality.CONVENTIONAL

        match = result.potential_matches[0]
        assert match.profile == "trp_6_7_FDH"
        assert match.confidence == 1.0  # TODO: but should it be?
        assert match.substrate == "tryptophan"
        assert match.target_positions == (7,)
        assert match.number_of_decorations == "mono"
        assert result.confidence == 3.

    def test_tryptophan_6(self):
        result = self.get_single_result("ktzR")

        assert result.confidence == 3.
        assert len(result.potential_matches) == 1
        assert result.conventionality == Conventionality.CONVENTIONAL

        match = result.potential_matches[0]
        assert match.profile == "trp_6_7_FDH"
        assert match.confidence == 1.0
        assert match.substrate == "tryptophan"
        assert match.target_positions == (6,)
        assert match.number_of_decorations == "mono"

    def test_tryptophan_5_weak(self):
        result = self.get_single_result("mibH")

        assert result.confidence == 2.5
        assert len(result.potential_matches) == 1
        assert result.conventionality == Conventionality.CONVENTIONAL

        match = result.potential_matches[0]
        assert match.profile == "trp_5_FDH"
        assert match.confidence == 0.5
        assert match.substrate == "tryptophan"
        assert match.target_positions == (5,)
        assert match.number_of_decorations == "mono"

    def test_tryptophan_5_strong(self):
        result = self.get_single_result("SpmH")

        assert result.confidence == 3.
        assert len(result.potential_matches) == 1
        assert result.conventionality == Conventionality.CONVENTIONAL

        match = result.potential_matches[0]
        assert match.profile == "trp_5_FDH"
        assert match.confidence == 1.
        assert match.substrate == "tryptophan"
        assert match.target_positions == (5,)
        assert match.number_of_decorations == "mono"

    def test_pyrrole_tetra(self):
        result = self.get_single_result("bmp2")

        assert result.confidence == 3.0
        assert len(result.potential_matches) == 1
        assert result.conventionality == Conventionality.CONVENTIONAL

        match = result.potential_matches[0]
        assert match.profile == "pyrrole_FDH"
        assert match.confidence == 1.0
        assert match.substrate == "pyrrole"
        assert match.target_positions == tuple()
        assert match.number_of_decorations == "tetra"

    def test_pyrrole_mono_di(self):
        result = self.get_single_result("HrmQ")

        assert result.confidence == 3.0
        assert len(result.potential_matches) == 1
        assert result.conventionality == Conventionality.CONVENTIONAL

        match = result.potential_matches[0]
        assert match.profile == "pyrrole_FDH"
        assert match.confidence == 1.0
        assert match.substrate == "pyrrole"
        assert match.target_positions == tuple()
        assert match.number_of_decorations == "mono_di"

    def test_orsellinic(self):
        result = self.get_single_result("XanH")

        assert result.confidence == 3.0
        assert result.conventionality == Conventionality.CONVENTIONAL

        assert len(result.potential_matches) == 1

        first = result.potential_matches[0]
        assert first.profile == "cycline_orsellinic_FDH"
        assert first.confidence == 1.0
        assert first.substrate == "cycline_orsellinic-like"
        assert first.target_positions == (6, 8)
        assert not first.number_of_decorations

    def test_hpg(self):
        result = self.get_single_result("End30")

        assert result.confidence == 3.0
        assert len(result.potential_matches) == 2
        assert result.conventionality == Conventionality.CONVENTIONAL
        assert result.cofactor == "flavin"
        assert result.family == "flavin-dependent"

        match = result.potential_matches[0]
        assert match.profile == "tyrosine-like_hpg_FDH"
        assert match.confidence == 1.0
        assert match.substrate == "Hpg"
        assert match.target_positions == (6, 8)
        assert not match.number_of_decorations

        assert result.potential_matches[1].substrate == "Tyr"
        assert match.confidence > result.potential_matches[1].confidence
