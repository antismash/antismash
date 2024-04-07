# License: GNU Affero General Public License v3 or later
# A copy of GNU AGPL v3 should have been included in this software package in LICENSE.txt.

# for test files, silence irrelevant and noisy pylint warnings
# pylint: disable=use-implicit-booleaness-not-comparison,protected-access,missing-docstring

from pathlib import Path
from typing import Union, Optional, List

from antismash.common.secmet import CDSFeature
from antismash.common.path import get_full_path
from antismash.common.signature import HmmSignature
from antismash.detection.genefunctions.halogenases.halogenases import (
    Match,
    HalogenaseHmmResult,
    TailoringEnzymes)
from antismash.detection.genefunctions.halogenases.flavin_dependent import substrate_analysis

SPECIFIC_PROFILES = [HmmSignature("pyrrole_FDH",
                                  "Pyrrole halogenase",
                                  400, get_full_path(str(Path(__file__).parents[1]),
                                                     "data", "pyrrole_FDH.hmm"))]

PYRROLE_SIGNATURE = [110, 111, 318, 322, 348, 362]

PYRROLE_SIGNATURE_RESIDUES = {"mono_di":"DRSVFW",
                              "unconv_mono_di":"YRRNFN",
                              "tetra":"RRYFFA"}

def search_for_match(name, residues, halogenase: TailoringEnzymes, hit: HalogenaseHmmResult,
                    position: Union[int, List[int]], cutoffs: List[float], *,
                    check_residues: bool = True, sig_residues: Union[str, dict[str,str]] = "",
                    confidence: float = 1, targets: bool = False):
    
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
    if name == "pyrrole_FDH":
        search_for_match(name, residues, halogenase, hit, 0,
                         cutoffs=[SPECIFIC_PROFILES[0].cutoff],
                         sig_residues=PYRROLE_SIGNATURE_RESIDUES,
                         targets = True)
        halogenase.substrates = "pyrrole"

def get_signatures() -> List[List[int]]:
    return [PYRROLE_SIGNATURE]

def get_consensus_signature(cds: CDSFeature, hit: HalogenaseHmmResult,
                 ) -> Union[dict, dict[str, str]]:
    
    if hit.query_id == "pyrrole_FDH":
        signature_residues: dict[str, Optional[str]] = {}
        substrates_signatures = dict(zip(list(PYRROLE_SIGNATURE_RESIDUES.keys()),(3*PYRROLE_SIGNATURE)))
        for substrate, signature in substrates_signatures.items():
            signature_residues[substrate] = substrate_analysis.search_residues(cds.translation,
                                                                               PYRROLE_SIGNATURE,
                                                                               hit)
        return {"pyrrole_FDH": signature_residues}
    return

