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

def search_for_match(name, residues, halogenase: TailoringEnzymes, hit: HalogenaseHmmResult,
                    position: Union[int, List[int]], cutoffs: List[float], *,
                    check_residues: bool = True, sig_residues: Union[str, dict[str,str]] = "",
                    confidence: float = 1):

    if hit.query_id != name:
        return False
    cutoffs.sort(reverse=True)
    modifier = 1.
    for cutoff in cutoffs:
        if hit.bitscore >= cutoff and (residues == sig_residues or not check_residues):
            halogenase.add_potential_matches(Match(hit.query_id, "flavin", "FDH", position,
                                                confidence * modifier, residues,""))
    return True

def update_match(name, residues, halogenase: TailoringEnzymes, hit: HalogenaseHmmResult) -> None:
    if name == "trp_5_FDH":
        search_for_match(name, residues, halogenase, hit, 5,
                        cutoffs=[SPECIFIC_PROFILES[0].cutoff,
                                 850],
                        sig_residues=TRP_5_SIGNATURE_RESIDUES)
        halogenase.substrates = "tryptophan"
    elif name == "trp_6_7_FDH":
        if not search_for_match(name, residues, halogenase, hit, 6, cutoffs=[770],
                               sig_residues=TRP_6_SIGNATURE_RESIDUES):
            search_for_match(name, residues, halogenase,
                            hit, 7, [SPECIFIC_PROFILES[1].cutoff], check_residues=False)
            halogenase.substrates = "tryptophan"

def get_signatures() -> List[List[int]]:
    return [TRP_5_SIGNATURE, TRP_6_SIGNATURE]

def get_consensus_signature(cds: CDSFeature, hit: HalogenaseHmmResult,
                 ) -> Union[dict, dict[str, str]]:
    residues = {}
    if hit.query_id == "trp_5_FDH":
        residues = substrate_analysis.retrieve_fdh_signature_residues(cds.translation, hit, TRP_5_SIGNATURE)
    if hit.query_id == "trp_6_7_FDH":
        residues = substrate_analysis.retrieve_fdh_signature_residues(cds.translation, hit, TRP_6_SIGNATURE)
    return residues
