# License: GNU Affero General Public License v3 or later
# A copy of GNU AGPL v3 should have been included in this software package in LICENSE.txt.

# for test files, silence irrelevant and noisy pylint warnings
# pylint: disable=use-implicit-booleaness-not-comparison,protected-access,missing-docstring

from pathlib import Path
from typing import Optional, List

from antismash.common.secmet import CDSFeature
from antismash.common.path import get_full_path
from antismash.common.signature import HmmSignature
from antismash.detection.genefunctions.halogenases.halogenases import (
    HalogenaseHmmResult,
    TailoringEnzymes)
from antismash.detection.genefunctions.halogenases.flavin_dependent import substrate_analysis

SPECIFIC_PROFILES = [HmmSignature("cycline_orsellinic_FDH",
                                  "Orsellinic acid-like or other phenolic substrate halogenase",
                                  500, get_full_path(Path(__file__).parents[2], "data",
                                                     "halogenases", "cycline_orsellinic_FDH.hmm")),
                     HmmSignature("tyrosine-like_hpg_FDH",
                                  "Tyrosine-like or Hpg substrate halogenase",
                                  300, get_full_path(Path(__file__).parents[2], "data",
                                                     "halogenases", "tyrosine-like_hpg_FDH.hmm"))]

TYROSINE_LIKE_SIGNATURE = [58, 74, 89, 92, 99, 107, 149, 150,
                            152, 209, 215, 217, 219, 245, 267,
                            268, 282, 284, 289, 290, 293, 295, 305, 331, 357]
HPG_SIGNATURE =  [66, 158, 196, 200, 246, 259]

OTHER_PHENOLIC_SIGNATURE = [23, 27, 39, 40, 59, 74, 109, 113, 120, 124, 133, 165,
                            166, 168, 233, 284, 291, 303, 305, 306, 309, 311]

TYROSINE_LIKE_SIGNATURE_RESIDUES = "GFQRLGDAGLSGVPSYGADPSGLYW"
HPG_SIGNATURE_RESIDUES = "SHCGMQ"

OTHER_PHENOLIC_SIGNATURE_RESIDUES = "LGPRGGRDAGVDAGGYGFDPSG"

def update_match(name, residues, halogenase: TailoringEnzymes, hit: HalogenaseHmmResult) -> None:
    if name == "tyrosine-like_hpg_FDH":
        if substrate_analysis.check_for_match(name, residues, halogenase, hit, [6, 8],
                            cutoffs=[300, 390],
                            sig_residues=TYROSINE_LIKE_SIGNATURE_RESIDUES):
            halogenase.substrate = "Tyr"
        if substrate_analysis.check_for_match(name, residues, halogenase, hit, [6, 8],
                            cutoffs=[300, 390],
                            sig_residues=HPG_SIGNATURE_RESIDUES):
            halogenase.substrate = "Hpg"
    elif name == "cycline_orsellinic_FDH":
        substrate_analysis.check_for_match(name, residues, halogenase, hit, [6, 8],
                                cutoffs=[500],
                                sig_residues=OTHER_PHENOLIC_SIGNATURE_RESIDUES)
        halogenase.substrate = "cycline_orsellinic"

def get_signatures() -> List[List[int]]:
    return [TYROSINE_LIKE_SIGNATURE,
            HPG_SIGNATURE,
            OTHER_PHENOLIC_SIGNATURE]

def get_consensus_signature(cds: CDSFeature, hit: HalogenaseHmmResult,
                 ) -> Optional[dict[str, Optional[str]]]:
    residues = {}
    if hit.query_id == "tyrosine-like_hpg_FDH":
        residues = substrate_analysis.get_residues(cds.translation, hit, [TYROSINE_LIKE_SIGNATURE,
                                                                          HPG_SIGNATURE],
                                 enzyme_substrates=["Tyr", "Hpg"])
        return {"tyrosine-like_hpg_FDH": residues}

    if hit.query_id == "cycline_orsellinic_FDH":
        residues = substrate_analysis.get_residues(cds.translation, hit,
                                                   OTHER_PHENOLIC_SIGNATURE)
    return residues
