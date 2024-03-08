# License: GNU Affero General Public License v3 or later
# A copy of GNU AGPL v3 should have been included in this software package in LICENSE.txt.

# for test files, silence irrelevant and noisy pylint warnings
# pylint: disable=use-implicit-booleaness-not-comparison,protected-access,missing-docstring

import logging
from typing import Any, Dict, List, Optional, Iterable, Union

from dataclasses import dataclass, field

from antismash.common.path import get_full_path
from antismash.common.hmmscan_refinement import HMMResult
from antismash.common.secmet import CDSFeature
from antismash.common.signature import HmmSignature
from antismash.common import (
    subprocessing,
    utils
)

NAME = "halogenases_analysis"
SHORT_DESCRIPTION = """Categorization of halogenases based on family and function.
                       It utilizes pHMMs with thresholds and the presence or absence of
                       characteristic/signature residues to identify the substrate or
                       the number of halogenation.
                       Signature residues for the pyrrole substrate were determined by
                       using information from an experimental study (https://doi.org/10.1073/pnas.1519695113).
                       Signature residues for all the other substrates were extracted computationally,
                       by looking at the positions that have the same value.
                    """

"""FX[DE]PX[]EFL"""
CHARGED_RESIDUES = [""]

GENERAL_FDH_PROFILES = [HmmSignature("all_general_FDH",
                                     "Member of the Flavin-dependent halogenase family",
                            100, get_full_path(__file__, "data",
                                               "halogenases", "all_general_FDH.hmm")),
                        HmmSignature("unconventional_FDH",
                                     "Unconventional member of the Flavin-dependent halogenase family",
                            100, get_full_path(__file__, "data",
                                               "halogenases", "unconventional_FDH.hmm"))]

SPECIFIC_FDH_PROFILES = [HmmSignature("trp_5_FDH",
                                      "Tryptophan-5 halogenase",
                                      350, get_full_path(__file__, "data",
                                                       "halogenases", "trp_5_FDH.hmm")),
                        HmmSignature("trp_6_7_FDH",
                                     "Tryptophan-6 or 7 halogenase",
                                     770, get_full_path(__file__, "data",
                                                       "halogenases", "trp_6_7_FDH.hmm")),
                        HmmSignature("pyrrole_FDH",
                                     "Pyrrole halogenase",
                                     400, get_full_path(__file__, "data",
                                                        "halogenases", "pyrrole_FDH.hmm")),
                        HmmSignature("cycline_orsellinic_FDH",
                                     "Orsellinic acid-like or other phenolic substrate halogenase",
                                     500, get_full_path(__file__, "data",
                                                        "halogenases", "cycline_orsellinic_FDH.hmm")),
                        HmmSignature("tyrosine-like_hpg_FDH",
                                     "Tyrosine-like or Hpg substrate halogenase",
                                     300, get_full_path(__file__, "data",
                                                        "halogenases", "tyrosine-like_hpg_FDH.hmm"))]

class HalogenaseHmmResult(HMMResult):
    """ Enzymes identified as a halogenase
        hit_id: name of the matching profile
        start: start position within the query's translation
        end: end position within the query's translation
        evalue: e-value of the hit
        bitscore: bitscore of the hit
        query_id: name of the profile
        enzyme_type: type of halogenase (e.g. Flavin-dependent, SAM-dependent)
        profile: path to the pHMM file
        internal_hits: any hits contained by this hit
    """
    def __init__(self, hit_id: str, bitscore: float, query_id: str, enzyme_type: str,
                 profile: str, start: int = 0, end: int = 0, evalue: float = 0.0,
                 internal_hits: Iterable[HMMResult] = None) -> None:
        super().__init__(hit_id, start, end, evalue, bitscore, internal_hits=internal_hits)
        self.query_id = query_id
        self.enzyme_type = enzyme_type
        self.profile = profile


@dataclass
class Match:
    """ Match of the enzyme categorized by check_for_fdh,
        with details about which pHMM (profile) was hit, what position
        the halogenation occurs, what is the confidence of the categorization,
        and what are the signature residues of the protein sequence"""
    profile: str
    position: Optional[Union[int, str, List[int]]]
    confidence: float
    signature: Optional[str]
    substrate: Optional[Union[int, str]] = None

    def to_json(self) -> dict[str, Any]:
        return vars(self)

    @classmethod
    def from_json(cls, data: dict[str, Any]) -> "Match":
        return cls(**data)


@dataclass
class HalogenaseResult:
    """ Details about the categorized enzymes
        cds_name: name of the cds that the result is for
        family: type of halogenase family (FD, SAMD, HD, VD, I/KGD)
        substrate: name of the common substrate (e.g. tryptophan, pyrrole, aliphatic molecule)
        position: position number of halogenated atom
        confidence: confidence of the categorization (between 0 and 1)
        signature: string of amino acid residues
        coenzyme: is coenzyme present or not
        potential_matches:  possible categorizations
                            (e.g. if enzyme meets requirements for more groups)"""
    cds_name: str
    family: str = "Halogenase"
    substrate: str = ""
    position: Optional[Union[int, str, List[int]]] = None
    confidence: float = 0
    signature: Optional[str] = None
    coenzyme: bool = False
    potential_matches: list[Match] = field(default_factory=list)

    def add_potential_matches(self, match: Match) -> None:
        """ Adds the features of an enzyme group to list"""
        self.potential_matches.append(match)

    @staticmethod
    def almost_equal(confidence: float, highest_confidence: float,
                     difference: float) -> bool:
        """ Helper function for get_best_match. Allows flexibility to
            handle confidence divergences."""
        return abs(confidence-highest_confidence) <= difference

    def get_best_match(self) -> list[Match]:
        """ If an enzyme meets the requirements for several groups,
            it compares the confidences of the categorizations and
            returns the one with the highest confidence or list of matches.
            If there are more groups with the same confidence, it returns the list of those."""
        best_match = []

        if self.potential_matches:
            if len(self.potential_matches) == 1:
                return [self.potential_matches[0]]

            highest_confidence = max(profile.confidence for profile in self.potential_matches)
            for profile in self.potential_matches:
                if self.almost_equal(profile.confidence, highest_confidence, 0.005):
                    best_match.append(profile)

        return best_match

    def finalize_enzyme(self) -> None:
        """ If there is a best match among the matches based on confidence,
            get that one match and define position, confidence and signature
            in the enzyme instance based on that.
            If there is no one best match, it doesn't change anything."""
        best_matches = self.get_best_match()
        assert isinstance(best_matches, list), best_matches
        if not best_matches:
            return

        if len(best_matches) == 1:
            best_match = best_matches[0]
            self.position = best_match.position
            self.confidence = best_match.confidence
            self.signature = best_match.signature
            if best_match.substrate:
                self.substrate = best_match.substrate

    def to_json(self) -> Dict[str, Any]:
        """ Constructs a JSON representation of this instance """
        potential_matches_json = [match.to_json() for match in self.potential_matches]

        return {
            "family": self.family,
            "gbk_id": self.cds_name,
            "substrate": self.substrate,
            "position": self.position,
            "confidence": self.confidence,
            "signature": self.signature,
            "coenzyme": self.coenzyme,
            "potential_matches": potential_matches_json
        }

    @classmethod
    def from_json(cls, data: Dict[str, Any]) -> "HalogenaseResult":
        """ Constructs the HalogenasesResults from the JSON representation """

        family = data["family"]
        gbk_id = data["gbk_id"]
        substrate = data["substrate"]
        position = data["position"]
        confidence = data["confidence"]
        signature = data["signature"]
        coenzyme = data["coenzyme"]
        matches = [Match.from_json(profile) for profile in data["potential_matches"]]
        enzyme = cls(gbk_id, family, substrate, position, confidence,
                     signature, coenzyme, matches)
        return enzyme


@dataclass
class HalogenaseCategories:
    name: str
    residues: Union[str, dict[str, str]]

    @staticmethod
    def get_signatures() -> List[List[int]]:
        raise NotImplementedError()

    def check_for_match(self, halogenase: HalogenaseResult, hit: HalogenaseHmmResult,
                        position: Union[int, List[int]], cutoffs: List[float], *,
                        check_residues: bool = True, sig_residues: Union[str, dict[str,str]] = "",
                        confidence: float = 1, targets: bool = False) -> bool:
        """in case the predefined signatures (self.signatures) is a dict,
            so there are several substrates, the searched signatures (sig_residues)
            also have to be a dict"""

        found = False
        if hit.query_id != self.name:
            return False
        cutoffs.sort(reverse=True)
        modifier = 1.
        for cutoff in cutoffs:
            if not isinstance(self.residues, dict):
                if hit.bitscore >= cutoff and (self.residues == sig_residues or not check_residues):
                    halogenase.add_potential_matches(Match(hit.query_id, position,
                                                        confidence * modifier, self.residues,""))
                    return True
                modifier = 0.5
            elif isinstance(self.residues, dict) and isinstance(sig_residues, dict):
                for subs, sig_res in self.residues.items():
                    if hit.bitscore >= cutoff and (sig_res == sig_residues[subs] or
                                                   not check_residues):
                        if not targets:
                            halogenase.add_potential_matches(Match(hit.query_id, position,
                                                                   confidence * modifier,
                                                                   sig_res, subs))
                        else:
                            halogenase.add_potential_matches(Match(hit.query_id, subs,
                                                                   confidence * modifier,
                                                                   sig_res))
                modifier = 0.5
                return True
            
            elif isinstance(self.residues, dict) and not isinstance(sig_residues, dict):
                for subs, sig_res in self.residues.items():
                    if hit.bitscore >= cutoff and (sig_res == sig_residues or
                                                   not check_residues):
                        if not targets:
                            halogenase.add_potential_matches(Match(hit.query_id, position,
                                                                   confidence * modifier,
                                                                   sig_res, subs))
                        else:
                            halogenase.add_potential_matches(Match(hit.query_id, subs,
                                                                   confidence * modifier,
                                                                   sig_res))
                modifier = 0.5
                return True                

        return found


class FlavinDependents(HalogenaseCategories):
    """ The second motif (WxWxI(P)) is responsible for preventing monooxygenating activity
        The third motif (Fx.Px.Sx.G) lines the tunnel between the flavin binding site and the site of halogenation"""
    def __init__(self, name: str, residues: dict[str, str]) -> None: # make conventional a property
        super().__init__(name, residues)
        self.family = "Flavin-dependent halogenases"
        self.profiles = ["general", "trp_5", "trp_6_7",
                         "pyrrole", "tyrosine-like_hpg", "cycline_orsellinic-like"]

    GENERAL_FDH_MOTIFS = {"W.W.I.": [206, 207, 208, 209, 210, 211],
                          "F.*P.*S.G": [280, 281, 282, 283, 284, 285, 286,
                                         287, 288, 289, 290, 291, 292, 293,
                                         294, 295, 296, 297, 298, 299, 300,
                                         301, 302, 303, 304]
                                         }

    TRP_5_SIGNATURE = [33, 35, 41, 75, 77, 78, 80, 92, 101, 128, 165, 186,
                       187, 195, 302, 310, 342, 350, 400, 446, 450, 452, 482]
    TRP_6_SIGNATURE = [19, 37, 45, 73, 75, 90, 129, 130, 142, 157, 181, 192, 194, 219,
                       221, 225, 227, 237, 287, 306, 337, 339, 350, 353, 356, 411, 462, 505]

    PYRROLE_SIGNATURE = [110, 111, 318, 322, 348, 362]

    TYROSINE_LIKE_SIGNATURE = [58, 74, 89, 92, 99, 107, 149, 150,
                               152, 209, 215, 217, 219, 245, 267,
                               268, 282, 284, 289, 290, 293, 295, 305, 331, 357]
    HPG_SIGNATURE =  [66, 158, 196, 200, 246, 259]

    OTHER_PHENOLIC_SIGNATURE = [23, 27, 39, 40, 59, 74, 109, 113, 120, 124, 133, 165,
                                166, 168, 233, 284, 291, 303, 305, 306, 309, 311]

    # signature residues
    TRP_5_SIGNATURE_RESIDUES = "VSILIREPGLPRGVPRAVLPGEA"
    TRP_6_SIGNATURE_RESIDUES = "TEGCAGFDAYHDRFGNADYGLSIIAKIL"

    PYRROLE_SIGNATURE_RESIDUES = {"mono_di":"DRSVFW",
                                  "unconv_mono_di":"YRRNFN",
                                  "tetra":"RRYFFA"}

    TYROSINE_LIKE_SIGNATURE_RESIDUES = "GFQRLGDAGLSGVPSYGADPSGLYW"
    HPG_SIGNATURE_RESIDUES = "SHCGMQ"

    OTHER_PHENOLIC_SIGNATURE_RESIDUES = "LGPRGGRDAGVDAGGYGFDPSG"

    @staticmethod
    def get_signatures() -> List[List[int]]:
        raise NotImplementedError()

class TryptophanSubstrate(FlavinDependents):
    def __init__(self, name: str, residues: str) -> None:
        super().__init__(name, residues)
        self.cutoff = None

    def update_match(self, halogenase: HalogenaseResult, hit: HalogenaseHmmResult) -> None:
        if self.name == "trp_5_FDH":
            self.check_for_match(halogenase, hit, 5,
                                 cutoffs=[350, 850],
                                 sig_residues=FlavinDependents.TRP_5_SIGNATURE_RESIDUES)
            halogenase.substrate = "tryptophan"
        elif self.name == "trp_6_7_FDH":
            if not self.check_for_match(halogenase, hit, 6, cutoffs=[770],
                                        sig_residues=FlavinDependents.TRP_6_SIGNATURE_RESIDUES):
                self.check_for_match(halogenase, hit, 7, [770], check_residues=False)
                halogenase.substrate = "tryptophan"

    @staticmethod
    def get_signatures() -> List[List[int]]:
        return [FlavinDependents.TRP_5_SIGNATURE, FlavinDependents.TRP_6_SIGNATURE]

    @staticmethod
    def get_residues(cds: CDSFeature, hit: HalogenaseHmmResult,
                     ) -> Optional[dict[str, Optional[str]]]:
        residues = None
        if hit.query_id == "trp_5_FDH":
            residues = get_residues(cds.translation, hit,
                                    FlavinDependents.TRP_5_SIGNATURE)
        if hit.query_id == "trp_6_7_FDH":
            residues = get_residues(cds.translation, hit,
                                    FlavinDependents.TRP_6_SIGNATURE)
        return residues

class PyrroleSubstrate(FlavinDependents):
    def __init__(self, name: str, residues: Union[str, dict]) -> None:
        super().__init__(name, residues)
        self.cutoff = 400
        self.signature_residues = FlavinDependents.PYRROLE_SIGNATURE_RESIDUES

    def update_match(self, halogenase: HalogenaseResult, hit: HalogenaseHmmResult) -> None:
        if self.name == "pyrrole_FDH":
            self.check_for_match(halogenase, hit, 0,
                                    cutoffs=[self.cutoff],
                                    sig_residues=self.signature_residues,
                                    targets = True)
            halogenase.substrate = "pyrrole"

    @staticmethod
    def get_signatures() -> List[List[int]]:
        return [FlavinDependents.PYRROLE_SIGNATURE]

    @staticmethod
    def get_residues(cds: CDSFeature, hit: HalogenaseHmmResult,
                     ) -> dict[str, Optional[str]]:
        residues = None
        if hit.query_id == "pyrrole_FDH":
            residues = get_residues(cds.translation, hit,
                                    PyrroleSubstrate.get_signatures(),
                                    substrates=["mono_di", "tetra", "unconv_mono_di"])
        return {"pyrrole_FDH": residues}

class PhenolicSubstrate(FlavinDependents):
    def __init__(self, name: str, residues: str) -> None:
        super().__init__(name, residues)
        self.cutoff = None
        self.other_phenolic = FlavinDependents.OTHER_PHENOLIC_SIGNATURE_RESIDUES

    def update_match(self, halogenase: HalogenaseResult, hit: HalogenaseHmmResult) -> None:
        if self.name == "tyrosine-like_hpg_FDH":
            if self.check_for_match(halogenase, hit, [6, 8],
                                cutoffs=[300, 390],
                                sig_residues=FlavinDependents.TYROSINE_LIKE_SIGNATURE_RESIDUES):
                halogenase.substrate = "Tyr"
            if self.check_for_match(halogenase, hit, [6, 8],
                                cutoffs=[300, 390],
                                sig_residues=FlavinDependents.HPG_SIGNATURE_RESIDUES):
                halogenase.substrate = "Hpg"
        elif self.name == "cycline_orsellinic_FDH":
            self.check_for_match(halogenase, hit, [6, 8],
                                    cutoffs=[500],
                                    sig_residues=self.other_phenolic)
            halogenase.substrate = "cycline_orsellinic"

    @staticmethod
    def get_signatures() -> List[List[int]]:
        return [FlavinDependents.TYROSINE_LIKE_SIGNATURE,
                FlavinDependents.HPG_SIGNATURE,
                FlavinDependents.OTHER_PHENOLIC_SIGNATURE]

    @staticmethod
    def get_residues(cds: CDSFeature, hit: HalogenaseHmmResult,
                     ) -> Optional[dict[str, Optional[str]]]:
        residues = {}
        if hit.query_id == "tyrosine-like_hpg_FDH":
            residues = get_residues(cds.translation, hit, [FlavinDependents.TYROSINE_LIKE_SIGNATURE,
                                     FlavinDependents.HPG_SIGNATURE],
                                     substrates=["Tyr", "Hpg"])
            return {"tyrosine-like_hpg_FDH": residues}

        elif hit.query_id == "cycline_orsellinic_FDH":
            residues = get_residues(cds.translation, hit,
                                    FlavinDependents.OTHER_PHENOLIC_SIGNATURE)
        return residues

def get_residues(translation: str, hmm_result: HalogenaseHmmResult,
                 signatures: Union[list[int], list[list[int]]], substrates: list = None)\
                 -> dict[str, Optional[str]]:
    """ Get signature residues for an enzyme from each pHMM

        Arguments:
            sequence: protein sequence
            hmm_result: instance of HmmResult class,
                        which contains information about the hit in a pHMM
            signatures: list of the positions that defines the signature residues in a pHMM

        Returns:
            signature residues which were retrieved from a certain pHMM
    """
    signature_residues: dict[str, Optional[str]] = {}

    if not substrates:
        residue = search_signature_residues(translation, signatures, hmm_result)
        signature_residues[hmm_result.query_id] = residue
    else:
        if len(signatures) == 1:
            substrates_signatures = dict(zip(substrates,(3*signatures)))
        else:
            substrates_signatures = dict(zip(substrates,signatures))
        for substrate, signature in substrates_signatures.items():
            signature_residues[substrate] = search_signature_residues(translation,
                                                                      signature, hmm_result)

    return signature_residues

def search_signature_residues(sequence: str, positions: list[int],
                              hmm_result: HalogenaseHmmResult,
                              max_evalue: float = 0.1) -> Optional[str]:
    """ Get the signature residues from the pHMM for the searched protein sequence

        Arguments:
            sequence: protein sequence
            positions: list of position numbers in the pHMM
            hmm_result: properties of the halogenase pHMM hit
            max_evalue: maximum e-value

        Returns:
            residues that are present in the given positions
    """
    args = ["-E", str(max_evalue)]

    results = subprocessing.hmmpfam.run_hmmpfam2(hmm_result.profile,
                                                 f">query\n{sequence}", extra_args=args)
    if not (results and results[0].hsps):
        logging.debug("no hits for query %s")
        return None

    found = False
    hit = None
    for hit in results[0].hsps:
        if hit.hit_id == hmm_result.query_id:
            found = True
            break

    if not found:
        logging.debug(
            "no hits for the enzyme in %s", hmm_result.query_id)
        return None

    profile = hit.aln[1].seq
    query = hit.aln[0].seq
    offset = hit.hit_start

    sites = utils.extract_by_reference_positions(query, profile,
                                                 [p - offset for p in positions if offset < p])

    return sites
    