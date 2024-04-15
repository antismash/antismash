# License: GNU Affero General Public License v3 or later
# A copy of GNU AGPL v3 should have been included in this software package in LICENSE.txt.

# for test files, silence irrelevant and noisy pylint warnings
# pylint: disable=use-implicit-booleaness-not-comparison,protected-access,missing-docstring

from pathlib import Path
from typing import Union, List

from antismash.common.secmet import CDSFeature
from antismash.common.path import get_full_path
from antismash.common.signature import HmmSignature
from antismash.detection.genefunctions.halogenases.halogenases import (
    Match,
    HalogenaseHmmResult,
    FlavinDependentHalogenases)
from antismash.detection.genefunctions.halogenases.flavin_dependent import substrate_analysis

SPECIFIC_PROFILES = [HmmSignature("tyrosine-like_hpg_FDH",
                                  "Tyrosine-like or Hpg substrate halogenase",
                                  390, get_full_path(str(Path(__file__).parents[1]),
                                                     "data", "tyrosine-like_hpg_FDH.hmm")),
                     HmmSignature("cycline_orsellinic_FDH",
                                  "Orsellinic acid-like or other phenolic substrate halogenase",
                                  500, get_full_path(str(Path(__file__).parents[1]),
                                                     "data", "cycline_orsellinic_FDH.hmm"))]

TYROSINE_LIKE_SIGNATURE = [58, 74, 89, 92, 99, 107, 149, 150,
                            152, 209, 215, 217, 219, 245, 267,
                            268, 282, 284, 289, 290, 293, 295, 305, 331, 357]
HPG_SIGNATURE =  [66, 158, 196, 200, 246, 259]

OTHER_PHENOLIC_SIGNATURE = [23, 27, 39, 40, 59, 74, 109, 113, 120, 124, 133, 165,
                            166, 168, 233, 284, 291, 303, 305, 306, 309, 311]

TYR_HPG_SIGNATURE_RESIDUES = {"Tyr": "GFQRLGDAGLSGVPSYGADPSGLYW",
                              "Hpg": "SHCGMQ"}

OTHER_PHENOLIC_SIGNATURE_RESIDUES = "LGPRGGRDAGVDAGGYGFDPSG"

def search_for_match(retrieved_residues: Union[dict[str, str], str], halogenase: FlavinDependentHalogenases,
                     hit: HalogenaseHmmResult, position: Union[int, List[int]],
                     cutoffs: Union[List[int], int], *, expected_residues: Union[str, dict[str, str]],
                     confidence: float = 1.) -> bool:
    """ Looks whether there are hmm hits that meet the requirement for the categorization

        Arguments:
            name: name of the substrate-specific pHMM
            residues: residues of the protein sequence in the place of the signature residues
            halogenase: initiated flavin-dependent halogenase
            hit: details of the hit (e.g. bitscore, name of the profile, etc.)
            position: position of decoration
            cutoffs: threshold(s) for the pHMM
            sig_residues: expected, substrate-specific signature residues
            confidence: reliability of the categorization

        Returns:
            if the hit is one of the tryptophan-specific pHMMs,
            then it adds the match, without returning anything,
            otherwise, it returns False
    """

    # check for halogenases with Tyr or Hpg substrates
    modifier = 1.
    if (isinstance(expected_residues, dict) and isinstance(cutoffs, list)
        and isinstance(retrieved_residues, dict)):
        cutoffs.sort(reverse=True)
        substrate_counter = 0
        for subs, sig_res in retrieved_residues.items():
            if sig_res == expected_residues[subs]:
                substrate_counter += 1

        for cutoff in cutoffs:
            # matches the residues for Tyrosine and Hpg as well
            if hit.bitscore < cutoff:
                modifier = .5
                continue
            if substrate_counter == 2:
                halogenase.add_potential_matches(Match(hit.query_id, "flavin", "FDH",
                                                    confidence * modifier, retrieved_residues["Hpg"],
                                                    target_positions=position, substrates="Hpg"))
                halogenase.add_potential_matches(Match(hit.query_id, "flavin", "FDH",
                                                       (confidence * modifier)-0.2, retrieved_residues["Tyr"],
                                                       target_positions=position, substrates="Tyr"))
                return True

            if retrieved_residues["Tyr"] == expected_residues["Tyr"]:
                halogenase.add_potential_matches(Match(hit.query_id, "flavin", "FDH",
                                                    confidence * modifier, retrieved_residues,
                                                    target_positions=position, substrates="Tyr"))
                return True
        return False
    if isinstance(cutoffs, int):
        if retrieved_residues != expected_residues or hit.bitscore < cutoffs:
            return False
        halogenase.add_potential_matches(Match(hit.query_id, "flavin", "FDH",
                                                    confidence * modifier, retrieved_residues,
                                                    target_positions=position,
                                                    substrates="cycline_orsellinic-like"))
        return True
    return False

def update_match(name: str, residues: dict[str, str], halogenase: FlavinDependentHalogenases,
                 hit: HalogenaseHmmResult) -> None:
    """ Looks whether there are hmm hits that meet the requirement for the categorization
        as Tyr, Hpg, or cycline/orsellinic-like halogenase

        Arguments:
            name: name of the substrate-specific pHMM
            residues: residues of the protein sequence in the place of the signature residues
            halogenase: initiated flavin-dependent halogenase
            hit: details of the hit (e.g. bitscore, name of the profile, etc.)

        Returns:
            if the categorization as Tyr/Hpg/other-phenolic-halogenase could be done
            it instanciates the match including the profile name, cofactor, family,
            position, confidence, signature and substrate,
            otherwise, it doesn't return anything and doesn't instanciate anything
    """

    if name == "tyrosine-like_hpg_FDH":
        search_for_match(residues, halogenase, hit, [6, 8],
                         cutoffs=[SPECIFIC_PROFILES[0].cutoff, 500],
                         expected_residues=TYR_HPG_SIGNATURE_RESIDUES)
    elif name == "cycline_orsellinic_FDH":
        search_for_match(residues, halogenase, hit, [6, 8],
                         cutoffs=SPECIFIC_PROFILES[1].cutoff,
                         expected_residues=OTHER_PHENOLIC_SIGNATURE_RESIDUES)

def get_consensus_signature(cds: CDSFeature, hit: HalogenaseHmmResult
                            ) -> Union[dict, dict[str, str]]:
    """ Retrieves the residues from the substrate-specific,
        pHMMs that are in the positions of the signature residues

        Arguments:
            cds: gene/CDS and its properties
            hit: details of the hit (e.g. bitscore, name of the profile, etc.)

        Returns:
            if the name of the pHMM doesn't match the substrate-specific one's,
            it returns an empty dictionary,
            otherwise, it returns the residues, that are in the same positions as
            the sugnature residues
    """

    residues = {}
    if hit.query_id == "tyrosine-like_hpg_FDH":
        residues = substrate_analysis.retrieve_fdh_signature_residues(cds.translation, hit,
                                                                      [TYROSINE_LIKE_SIGNATURE,
                                                                       HPG_SIGNATURE],
                                                                      enzyme_substrates=["Tyr", "Hpg"])
        return {"tyrosine-like_hpg_FDH": residues}

    if hit.query_id == "cycline_orsellinic_FDH":
        residues = substrate_analysis.retrieve_fdh_signature_residues(cds.translation, hit,
                                                                      OTHER_PHENOLIC_SIGNATURE)
    return residues
