# License: GNU Affero General Public License v3 or later
# A copy of GNU AGPL v3 should have been included in this software package in LICENSE.txt.

# for test files, silence irrelevant and noisy pylint warnings
# pylint: disable=use-implicit-booleaness-not-comparison,protected-access,missing-docstring

import os
import unittest
from unittest.mock import patch

from antismash.detection.genefunctions.tools.halogenases.data_structures import (
    HalogenaseHit as HalogenaseHmmResult,
    Match,
    Profile,
)
from antismash.detection.genefunctions.tools.halogenases.flavin_dependent import (
    substrate_analysis,
)
from antismash.detection.genefunctions.tools.halogenases.flavin_dependent.substrate_analysis import (
    SUBGROUPS,
    categorise_on_substrate_level,
)
FDH = HalogenaseHmmResult

indolic = SUBGROUPS["indolic"]

TRP_5 = indolic.profiles["Trp_5"]
TRP_6_7 = indolic.profiles["Trp_6_7"]
TRP_6_7_PROFILE_NAME = TRP_6_7.profile_name
TRP_6, TRP_7 = TRP_6_7.motifs


def create_hit(**kwargs):
    options = {
        "query_id": "query",
        "reference_id": "reference",
        "bitscore": 50,
        "evalue": 1e-20,
        "conventionality_residues": "ABCDEF",
        "profile": Profile(name="name", description="desc", profile_name=kwargs.get("reference_id", "reference"),
                           cofactor="flavin", family="flavin-dependent", filename=os.devnull, cutoffs=(5, 10)),
    }
    options.update(kwargs)
    return FDH(**options)


def create_flavin_match(profile, confidence=0., consensus_residues="", substrate=None,
                        target_positions=None, number_of_decorations=None,
                        ):
    return Match(
        profile,
        "flavin",
        "flavin-dependent",
        confidence=confidence,
        consensus_residues=consensus_residues,
        substrate=substrate,
        target_positions=target_positions,
        number_of_decorations=number_of_decorations,
    )


def create_motif_residue_mapping(profile):
    return {motif.name: motif.residues for motif in profile.motifs}


class IndolicBase(unittest.TestCase):
    def setUp(self):
        self.test_trp_5_match = create_flavin_match(
            TRP_5, confidence=1, substrate="tryptophan", target_positions=5,
            number_of_decorations="mono",
        )
        self.test_trp_6_7_match = create_flavin_match(
            TRP_6_7_PROFILE_NAME, confidence=1, substrate="tryptophan",
            target_positions=6, number_of_decorations="mono",
        )

        self.trp_5_hmm_result = HalogenaseHmmResult(
            query_id="query",
            bitscore=1000,
            evalue=1e-20,
            reference_id=TRP_5.profile_name,
            profile=TRP_5.profile_name,
        )
        self.trp_6_7_hmm_result = HalogenaseHmmResult(
            query_id="query",
            bitscore=1000,
            evalue=1e-20,
            reference_id=TRP_6_7_PROFILE_NAME,
            profile=TRP_6_7_PROFILE_NAME,
        )

        tryptophan_single_matches = [self.test_trp_6_7_match]
        tryptophan_matches = [self.test_trp_5_match, self.test_trp_6_7_match]

        # Trp-5 halogenase
        self.trp_5_enzyme_with_matches = create_hit(query_id="mibH", potential_matches=tryptophan_matches)
        # Trp-6 halogenase
        self.trp_enzyme_with_matches = create_hit(query_id="ktzR", potential_matches=tryptophan_single_matches)
        # Trp-7 halogenase
        self.trp_with_no_matches = create_hit()


class TestIndolic(IndolicBase):
    def test_strong_trp_5(self):
        with patch.object(substrate_analysis, "retrieve_fdh_signature_residues",
                          return_value=create_motif_residue_mapping(TRP_5)) as patched:
            matches = categorise_on_substrate_level("MAGIC", [self.trp_5_hmm_result])
            patched.assert_called_once_with(
                "MAGIC", self.trp_5_hmm_result, TRP_5.motifs,
            )

        match = matches[0]
        assert match.profile == TRP_5.profile_name
        assert match.confidence == 1
        assert match.substrate == "tryptophan"
        assert match.number_of_decorations == "mono"
        assert match.target_positions == (5,)

    def test_weak_trp_5(self):
        low_quality_hit = HalogenaseHmmResult(reference_id=TRP_5.profile_name, bitscore=380,
                                              evalue=1e-2, query_id="test", profile=TRP_5.profile_name)
        with patch.object(substrate_analysis, "retrieve_fdh_signature_residues",
                          return_value=create_motif_residue_mapping(TRP_5)):
            matches = categorise_on_substrate_level("MAGIC", [low_quality_hit])
        match = matches[0]
        assert match.profile == TRP_5.profile_name
        assert match.confidence == 0.5
        assert match.substrate == "tryptophan"
        assert match.number_of_decorations == "mono"
        assert match.target_positions == (5,)

    def test_trp_6(self):
        self.trp_6_7_hmm_result.profile = TRP_6_7
        with patch.object(substrate_analysis, "retrieve_fdh_signature_residues",
                          return_value=create_motif_residue_mapping(TRP_6_7)):
            matches = categorise_on_substrate_level("MAGIC", [self.trp_6_7_hmm_result])
        match = matches[0]
        assert match.profile == TRP_6_7_PROFILE_NAME
        assert match.confidence == 1.0
        assert match.substrate == "tryptophan"
        assert match.number_of_decorations == "mono"
        assert match.target_positions == (6,)

    def test_trp_7(self):
        self.trp_6_7_hmm_result.profile = TRP_7
        residues = create_motif_residue_mapping(TRP_6_7)
        # add at least 5 different positions
        distant = ("-" * TRP_7.minimum_distance) + residues[TRP_7.name][TRP_7.minimum_distance:]
        for key in residues:
            residues[key] = distant
        with patch.object(substrate_analysis, "retrieve_fdh_signature_residues",
                          return_value=residues):
            matches = categorise_on_substrate_level("MAGIC", [self.trp_6_7_hmm_result])
        match = matches[0]
        assert match.profile == TRP_6_7_PROFILE_NAME
        assert match.confidence == 1.
        assert match.target_positions == (7,)
        assert match.substrate == "tryptophan"
        assert match.number_of_decorations == "mono"

    @patch.object(substrate_analysis, "retrieve_fdh_signature_residues",
                  return_value=create_motif_residue_mapping(TRP_5))
    def test_categorise_substrate_no_match(self, _patched_retrieve_fdh_signature_residues):
        result = categorise_on_substrate_level("MAGIC", [])
        assert not result

    @patch.object(substrate_analysis, "retrieve_fdh_signature_residues",
                  return_value=create_motif_residue_mapping(TRP_5))
    def test_categorise_substrate_good_match(self, _patched_consensus_sig):
        matches = categorise_on_substrate_level("MAGIC", [self.trp_5_hmm_result])
        assert matches == [
            create_flavin_match(profile=TRP_5.profile_name, confidence=1.0,
                                consensus_residues=TRP_5.motifs[0].residues,
                                substrate="tryptophan", target_positions=(5,),
                                number_of_decorations="mono")
        ]
