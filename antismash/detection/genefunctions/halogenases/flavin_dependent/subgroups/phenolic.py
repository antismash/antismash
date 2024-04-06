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
    TailoringEnzymes)
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

TYROSINE_LIKE_SIGNATURE_RESIDUES = "GFQRLGDAGLSGVPSYGADPSGLYW"
HPG_SIGNATURE_RESIDUES = "SHCGMQ"

OTHER_PHENOLIC_SIGNATURE_RESIDUES = "LGPRGGRDAGVDAGGYGFDPSG"

def search_for_match(name, residues, halogenase: TailoringEnzymes, hit: HalogenaseHmmResult,
                    position: Union[int, List[int]], cutoffs: List[float], *,
                    check_residues: bool = True, sig_residues: Union[str, dict[str,str]] = "",
                    confidence: float = 1):

    if hit.query_id != name:
        return False
    cutoffs.sort(reverse=True)
    modifier = 1.
    for cutoff in cutoffs:
        if isinstance(residues, dict) and isinstance(sig_residues, dict):
            for subs, sig_res in residues.items():
                if hit.bitscore >= cutoff and (sig_res == sig_residues or
                                                not check_residues):
                        halogenase.add_potential_matches(Match(hit.query_id, "flavin", "FDH",
                                                                position, confidence * modifier,
                                                                sig_res, subs))
                else:
                    halogenase.add_potential_matches(Match(hit.query_id, "flavin", "FDH", subs,
                                                            confidence * modifier,
                                                            sig_res))
            return True
        if hit.bitscore >= cutoff and (residues == sig_residues or not check_residues):
                halogenase.add_potential_matches(Match(hit.query_id, "flavin", "FDH", position,
                                                    confidence * modifier, residues,""))
                return True

def update_match(name, residues, halogenase: TailoringEnzymes, hit: HalogenaseHmmResult) -> None:
    if name == "tyrosine-like_hpg_FDH":
        if search_for_match(name, residues, halogenase, hit, [6, 8],
                            cutoffs=[SPECIFIC_PROFILES[0].cutoff, 500],
                            sig_residues=TYROSINE_LIKE_SIGNATURE_RESIDUES):
            halogenase.substrates = "Tyr"
        if search_for_match(name, residues, halogenase, hit, [6, 8],
                            cutoffs=[SPECIFIC_PROFILES[0].cutoff, 500],
                            sig_residues=HPG_SIGNATURE_RESIDUES):
            halogenase.substrates = "Hpg"
    elif name == "cycline_orsellinic_FDH":
        search_for_match(name, residues, halogenase, hit, [6, 8],
                         cutoffs=[SPECIFIC_PROFILES[1].cutoff],
                         sig_residues=OTHER_PHENOLIC_SIGNATURE_RESIDUES)
        halogenase.substrates = "cycline_orsellinic"

def get_signatures() -> List[List[int]]:
    return [TYROSINE_LIKE_SIGNATURE,
            HPG_SIGNATURE,
            OTHER_PHENOLIC_SIGNATURE]

def get_consensus_signature(cds: CDSFeature, hit: HalogenaseHmmResult,
                 ) -> Union[dict, dict[str, str]]:
    residues = {}
    if hit.query_id == "tyrosine-like_hpg_FDH":
        residues = substrate_analysis.retrieve_signature_residues(cds.translation, hit,
                                                   [TYROSINE_LIKE_SIGNATURE,
                                                    HPG_SIGNATURE],
                                                   enzyme_substrates=["Tyr", "Hpg"])
        return {"tyrosine-like_hpg_FDH": residues}

    if hit.query_id == "cycline_orsellinic_FDH":
        residues = substrate_analysis.retrieve_signature_residues(cds.translation, hit,
                                                   OTHER_PHENOLIC_SIGNATURE)
    return residues
