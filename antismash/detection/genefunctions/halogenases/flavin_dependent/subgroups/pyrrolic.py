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
    FlavinDependentHalogenases)
from antismash.detection.genefunctions.halogenases.flavin_dependent import substrate_analysis

SPECIFIC_PROFILES = [HmmSignature("pyrrole_FDH",
                                  "Pyrrole halogenase",
                                  400, get_full_path(str(Path(__file__).parents[1]),
                                                     "data", "pyrrole_FDH.hmm"))]

PYRROLE_SIGNATURE = [110, 111, 318, 322, 348, 362]

PYRROLE_SIGNATURE_RESIDUES = {"mono_di":"DRSVFW",
                              "unconv_mono_di":"YRRNFN",
                              "tetra":"RRYFFA"}

def search_for_match(residues: dict[str, str], halogenase: FlavinDependentHalogenases,
                     hit: HalogenaseHmmResult, cutoff: float, *,
                     sig_residues: Union[str, dict[str,str]] = "", confidence: float = 1
                     ) -> bool:
    """ Looks whether there are hmm hits that meet the requirement for the categorization

        Arguments:
            name: name of the substrate-specific pHMM
            residues: residues of the protein sequence in the place of the signature residues
            halogenase: initiated flavin-dependent halogenase
            hit: details of the hit (e.g. bitscore, name of the profile, etc.)
            cutoffs: threshold(s) for the pHMM
            sig_residues: substrate-specific signature residues
            confidence: reliability of the categorization

        Returns:
            if the hit is one of the tryptophan-specific pHMMs,
            then it adds the match, without returning anything,
            otherwise, it returns nothing
    """
    if (not isinstance(sig_residues, dict) or hit.bitscore < cutoff):
        return False
    for subs, sig_res in sig_residues.items():
        if residues == sig_residues[subs]:
            halogenase.add_potential_matches(Match(hit.query_id, "flavin", "FDH",
                                                    confidence, sig_res,
                                                    number_of_decorations=subs,
                                                    substrates = "pyrrole"))
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
    if name == "pyrrole_FDH":
        search_for_match(residues, halogenase, hit,
                         cutoff=SPECIFIC_PROFILES[0].cutoff,
                         sig_residues=PYRROLE_SIGNATURE_RESIDUES)

def get_consensus_signature(cds: CDSFeature, hit: HalogenaseHmmResult
                            ) -> dict[str, Optional[str]]:
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

    if hit.query_id == "pyrrole_FDH":
        signature_residues = substrate_analysis.search_residues(cds.translation,
                                                                PYRROLE_SIGNATURE,
                                                                hit)
    return {"pyrrole_FDH": signature_residues}
