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

SPECIFIC_PROFILES = [HmmSignature("trp_5_FDH",
                                 "Tryptophan-5 halogenase",
                                 350, get_full_path(Path(__file__).parents[2], "data",
                                                    "halogenases", "trp_5_FDH.hmm")),
                    HmmSignature("trp_6_7_FDH",
                                 "Tryptophan-6 or 7 halogenase",
                                 770, get_full_path(Path(__file__).parents[2], "data",
                                                    "halogenases", "trp_6_7_FDH.hmm"))]

TRP_5_SIGNATURE = [33, 35, 41, 75, 77, 78, 80, 92, 101, 128, 165, 186,
                    187, 195, 302, 310, 342, 350, 400, 446, 450, 452, 482]
TRP_6_SIGNATURE = [19, 37, 45, 73, 75, 90, 129, 130, 142, 157, 181, 192, 194, 219,
                    221, 225, 227, 237, 287, 306, 337, 339, 350, 353, 356, 411, 462, 505]

TRP_5_SIGNATURE_RESIDUES = "VSILIREPGLPRGVPRAVLPGEA"
TRP_6_SIGNATURE_RESIDUES = "TEGCAGFDAYHDRFGNADYGLSIIAKIL"

def update_match(name, residues, halogenase: TailoringEnzymes, hit: HalogenaseHmmResult) -> None:
    if name == "trp_5_FDH":
        substrate_analysis.check_for_match(name, residues, halogenase, hit, 5,
                                cutoffs=[350, 850],
                                sig_residues=TRP_5_SIGNATURE_RESIDUES)
        halogenase.substrates = "tryptophan"
    elif name == "trp_6_7_FDH":
        if not substrate_analysis.check_for_match(name, residues, halogenase, hit, 6, cutoffs=[770],
                                    sig_residues=TRP_6_SIGNATURE_RESIDUES):
            substrate_analysis.check_for_match(name, residues, halogenase,
                                               hit, 7, [770], check_residues=False)
            halogenase.substrates = "tryptophan"

def get_signatures() -> List[List[int]]:
    return [TRP_5_SIGNATURE, TRP_6_SIGNATURE]

def get_consensus_signature(cds: CDSFeature, hit: HalogenaseHmmResult,
                 ) -> Optional[dict[str, Optional[str]]]:
    residues = None
    if hit.query_id == "trp_5_FDH":
        residues = substrate_analysis.get_residues(cds.translation, hit, TRP_5_SIGNATURE)
    if hit.query_id == "trp_6_7_FDH":
        residues = substrate_analysis.get_residues(cds.translation, hit, TRP_6_SIGNATURE)
    return residues
