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

SPECIFIC_PROFILES = [HmmSignature("trp_5_FDH",
                                 "Tryptophan-5 halogenase",
                                 350, get_full_path(str(Path(__file__).parents[1]),
                                                    "data", "trp_5_FDH.hmm")),
                    HmmSignature("trp_6_7_FDH",
                                 "Tryptophan-6 or 7 halogenase",
                                 770, get_full_path(str(Path(__file__).parents[1]),
                                                    "data", "trp_6_7_FDH.hmm"))]

TRP_5_SIGNATURE = [33, 35, 41, 75, 77, 78, 80, 92, 101, 128, 165, 186,
                    187, 195, 302, 310, 342, 350, 400, 446, 450, 452, 482]
TRP_6_SIGNATURE = [19, 37, 45, 73, 75, 90, 129, 130, 142, 157, 181, 192, 194, 219,
                    221, 225, 227, 237, 287, 306, 337, 339, 350, 353, 356, 411, 462, 505]

TRP_5_SIGNATURE_RESIDUES = "VSILIREPGLPRGVPRAVLPGEA"
TRP_6_SIGNATURE_RESIDUES = "TEGCAGFDAYHDRFGNADYGLSIIAKIL"

def search_for_match(retrieved_residues: str, halogenase: FlavinDependentHalogenases,
                     hit: HalogenaseHmmResult, position: Union[int, List[int]],
                     cutoffs: List[float], *, check_residues: bool = True,
                     expected_residues: Union[str, dict[str,str]] = "", confidence: float = 1) -> bool:
    """ Looks whether there are hmm hits that meet the requirement for the categorization

        Arguments:
            retrieved_residues: residues of the protein sequence in the place of the signature residues
            halogenase: initiated flavin-dependent halogenase
            hit: details of the hit (e.g. bitscore, name of the profile, etc.)
            position: position of decoration
            cutoffs: threshold(s) for the pHMM
            check_residues: should the signature residues be looked at or not
            expected_residues: substrate-specific signature residues
            confidence: reliability of the categorization

        Returns:
            if the hit is one of the tryptophan-specific pHMMs,
            then it adds the match, without returning anything,
            otherwise, it returns False
    """
    cutoffs.sort(reverse=True)
    modifier = 1.
    for cutoff in cutoffs:
        if hit.bitscore < cutoff:
            modifier = .5
            continue
        if retrieved_residues == expected_residues or not check_residues:
            halogenase.add_potential_matches(Match(hit.query_id,"flavin", "FDH",
                                                   confidence * modifier, retrieved_residues,
                                                   target_positions=position,
                                                   number_of_decorations="mono",
                                                   substrates="tryptophan"))
            return True
    return False

def update_match(name: str, retrieved_residues: str, halogenase: FlavinDependentHalogenases,
                 hit: HalogenaseHmmResult) -> None:
    """ Looks whether there are hmm hits that meet the requirement for the categorization
        as Trp-5, Trp-6, or Trp-7 halogenase

        Arguments:
            name: name of the substrate-specific pHMM
            retrieved_residues: residues of the protein sequence in the place of the signature residues
            halogenase: initiated flavin-dependent halogenase
            hit: details of the hit (e.g. bitscore, name of the profile, etc.)

        Returns:
            if the categorization as Trp-5/6/7-halogenase could be done, it instanciates the match
            including the profile name, cofactor, family, position,
            confidence, signature and substrate,
            otherwise, it doesn't return anything and doesn't instanciate anything
    """
    if name == "trp_5_FDH":
        if search_for_match(retrieved_residues, halogenase, hit, 5,
                            cutoffs=[SPECIFIC_PROFILES[0].cutoff, 850],
                            expected_residues=TRP_5_SIGNATURE_RESIDUES):
            return
    elif name == "trp_6_7_FDH":
        if not search_for_match(retrieved_residues, halogenase, hit, 6, cutoffs=[770],
                               expected_residues=TRP_6_SIGNATURE_RESIDUES):
            if search_for_match(retrieved_residues, halogenase,
                                hit, 7, [SPECIFIC_PROFILES[1].cutoff], check_residues=False):
                return

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
            the signature residues
    """
    residues = {}
    if hit.query_id == "trp_5_FDH":
        residues = substrate_analysis.retrieve_fdh_signature_residues(cds.translation,
                                                                      hit, TRP_5_SIGNATURE)
    if hit.query_id == "trp_6_7_FDH":
        residues = substrate_analysis.retrieve_fdh_signature_residues(cds.translation,
                                                                      hit, TRP_6_SIGNATURE)
    return residues
