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

SPECIFIC_PROFILES = [HmmSignature("pyrrole_FDH",
                                  "Pyrrole halogenase",
                                  400, get_full_path(Path(__file__).parents[2], "data",
                                                     "halogenases", "pyrrole_FDH.hmm"))]

PYRROLE_SIGNATURE = [110, 111, 318, 322, 348, 362]

PYRROLE_SIGNATURE_RESIDUES = {"mono_di":"DRSVFW",
                              "unconv_mono_di":"YRRNFN",
                              "tetra":"RRYFFA"}

CUTOFF = 400

def update_match(name, residues, halogenase: TailoringEnzymes, hit: HalogenaseHmmResult) -> None:
    if name == "pyrrole_FDH":
        substrate_analysis.check_for_match(name, residues, halogenase, hit, 0,
                                cutoffs=[CUTOFF],
                                sig_residues=PYRROLE_SIGNATURE_RESIDUES,
                                targets = True)
        halogenase.substrates = "pyrrole"

def get_signatures() -> List[List[int]]:
    return [PYRROLE_SIGNATURE]

def get_consensus_signature(cds: CDSFeature, hit: HalogenaseHmmResult,
                 ) -> dict[str, Optional[str]]:
    residues = None
    if hit.query_id == "pyrrole_FDH":
        residues = substrate_analysis.get_residues(cds.translation, hit,
                                get_signatures(),
                                enzyme_substrates=["mono_di", "tetra", "unconv_mono_di"])
    return {"pyrrole_FDH": residues}
