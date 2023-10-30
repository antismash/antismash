# License: GNU Affero General Public License v3 or later
# A copy of GNU AGPL v3 should have been included in this software package in LICENSE.txt.

# for test files, silence irrelevant and noisy pylint warnings
# pylint: disable=use-implicit-booleaness-not-comparison,protected-access,missing-docstring

import unittest
from unittest.mock import Mock, patch

from antismash.common import fasta, path, subprocessing
from antismash.common.test.helpers import (
    FakeHSPHit,
    FakeHit,
)
from antismash.detection.genefunctions.tools.halogenases.data_structures import (
    HalogenaseHit,
)
from antismash.detection.genefunctions.tools.halogenases.flavin_dependent import (
    substrate_analysis,
)
from antismash.detection.genefunctions.tools.halogenases.flavin_dependent.substrate_analysis import (
    SUBGROUPS,
    categorise_on_substrate_level,
)
FDH = HalogenaseHit

phenolic = SUBGROUPS["phenolic"]
TYR_HPG = phenolic.profiles["tyrosine-hpg"]
ORSELLINIC = phenolic.profiles["orsellinic"]

TRANSLATIONS = fasta.read_fasta(path.get_full_path(__file__, "data", "translations.fasta"))
TRANSLATIONS = {key.rsplit("|", 1)[-1]: value for key, value in TRANSLATIONS.items()}


def create_motif_residue_mapping(profile):
    return {motif.name: motif.residues for motif in profile.motifs}


class TestPhenolic(unittest.TestCase):
    def setUp(self):
        self.hpg_hmm_result = HalogenaseHit(
            reference_id=TYR_HPG.profile_name,
            query_id="query_name",
            bitscore=600,
            evalue=1e-20,
            profile=TYR_HPG,
        )
        self.cycline_orsellinic_hmm_result = HalogenaseHit(
            reference_id=ORSELLINIC.profile_name,
            query_id="query_name",
            bitscore=600,
            evalue=1e-20,
            profile=ORSELLINIC,
        )

    def test_categorising_orsellinic(self):
        with patch.object(substrate_analysis, "retrieve_fdh_signature_residues",
                          return_value=create_motif_residue_mapping(ORSELLINIC)):
            matches = categorise_on_substrate_level("MAGIC", [self.cycline_orsellinic_hmm_result])
        assert len(matches) == 1
        match = matches[0]
        assert match.profile == ORSELLINIC.profile_name
        assert match.confidence == 1
        assert match.substrate == ORSELLINIC.motifs[0].substrate
        assert match.target_positions == (6, 8)

    def test_no_result(self):
        with patch.object(substrate_analysis, "extract_residues",
                          return_value=""):
            result = categorise_on_substrate_level("MAGIC", [self.cycline_orsellinic_hmm_result])
        assert not result

    def test_both_tyr_and_hpg(self):
        with patch.object(substrate_analysis, "retrieve_fdh_signature_residues",
                          return_value=create_motif_residue_mapping(TYR_HPG)):
            matches = categorise_on_substrate_level("MAGIC", [self.hpg_hmm_result])

        hpg_match, tyr_match = sorted(matches, key=lambda x: x.substrate)

        assert hpg_match.profile == TYR_HPG.profile_name
        assert hpg_match.confidence == 1.
        assert hpg_match.substrate == "Hpg"
        assert hpg_match.target_positions == (6, 8)

        assert tyr_match.profile == TYR_HPG.profile_name
        assert tyr_match.confidence == 0.8
        assert tyr_match.substrate == "Tyr"
        assert tyr_match.target_positions == (6, 8)

    def test_weak_hit_good_motif(self):
        hit = HalogenaseHit(
            query_id="query",
            reference_id=TYR_HPG.profile_name,
            bitscore=310,  # weak score
            evalue=1e-10,
            profile="dummy_path",
        )
        motifs_tyrosine_match_not_hpg = {
            "Tyr": "GFQRLGDAGLSGVPSYGADPSGLYW",
            "Hpg": "VALAMI",
        }
        with patch.object(substrate_analysis, "retrieve_fdh_signature_residues",
                          return_value=motifs_tyrosine_match_not_hpg):
            matches = categorise_on_substrate_level("MAGIC", [hit])
        match = matches[0]
        assert match.profile == TYR_HPG.profile_name
        assert match.confidence == 0.5
        assert match.substrate == "Tyr"
        assert match.target_positions == (6, 8)

    def test_good_hit_bad_motif(self):
        hit = HalogenaseHit(
            query_id="query",
            reference_id=TYR_HPG.profile_name,
            bitscore=1000,
            evalue=1e-20,
            profile=TYR_HPG.filename,
        )
        with patch.object(substrate_analysis, "retrieve_fdh_signature_residues",
                          return_value={}):
            matches = categorise_on_substrate_level("MAGIC", [hit])
        assert not matches

    def test_retrieve_fdh_signature_residues(self):
        alignment = FakeHSPHit(hit_id=TYR_HPG.profile_name, query_id="query_name")
        alignment.aln = [Mock(seq=TRANSLATIONS["BhaA"]), Mock(seq=TRANSLATIONS["BhaA"])]
        hit = FakeHit(1, 2, 1000, "foo")
        hit.hsps = [alignment]

        with patch.object(subprocessing.hmmpfam, "run_hmmpfam2", return_value=[hit]):
            residues = substrate_analysis.retrieve_fdh_signature_residues(
                "MAGIC", self.hpg_hmm_result, TYR_HPG.motifs,
            )
        assert residues == {
            "Tyr": "MRFGGVTDDNFSWQADVKQQYRANV",
            "Hpg": "AVFKGC",
        }
