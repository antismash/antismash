# License: GNU Affero General Public License v3 or later
# A copy of GNU AGPL v3 should have been included in this software package in LICENSE.txt.

# for test files, silence irrelevant and noisy pylint warnings
# pylint: disable=use-implicit-booleaness-not-comparison,protected-access,missing-docstring

import logging
from collections import defaultdict
from typing import Any, Dict, List, Optional, Iterable

from dataclasses import dataclass, field

from antismash.common.path import get_full_path
from antismash.common.hmmscan_refinement import HMMResult
from antismash.common.secmet import Record, CDSFeature
from antismash.common.signature import HmmSignature
from antismash.common import (
    subprocessing,
    path,
    fasta,
    utils,
)

NAME = "halogenases_analysis"
SHORT_DESCRIPTION = """categorization of halogenases based on family and function"""

pHMM_SIGNATURES = [
    HmmSignature("trp_5", "Flavin-dependent halogenase, modifying the 5th position of a Trp",
                 350, get_full_path(__file__, "data", "halogenases", "trp_5_v2.hmm")),
    HmmSignature("trp_6_7", "Flavin-dependent halogenase, modifying the 6/7th position of a Trp",
                 770, get_full_path(__file__, "data", "halogenases", "trp_6_7_v2.hmm"))]


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
    position: int
    confidence: float
    signature: Optional[str]

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
    position: Optional[int] = None
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
    residues: str

    @staticmethod
    def get_signatures() -> List[List[int]]:
        raise NotImplementedError()

    def check_for_match(self, halogenase: HalogenaseResult, hit: HalogenaseHmmResult,
                        position: int, cutoffs: List[float], *,
                        check_residues: bool = True, sig_residues: Optional[str] = None,
                        confidence: float = 1) -> bool:
        found = False
        if hit.query_id != self.name:
            return False

        cutoffs.sort(reverse=True)
        modifier = 1.
        for cutoff in cutoffs:
            if hit.bitscore >= cutoff and (self.residues == sig_residues or not check_residues):
                halogenase.add_potential_matches(Match(hit.hit_id, position,
                                                       confidence * modifier, self.residues))
                return True
            modifier = 0.5

        return found


class FlavinDependents(HalogenaseCategories):
    def __init__(self, name: str, residues: str) -> None:
        super().__init__(name, residues)
        self.family = "Flavin-dependent halogenases"

    TRP_5_SIGNATURE = [33, 35, 41, 75, 77, 78, 80, 92, 101, 128, 165, 186,
                       187, 195, 302, 310, 342, 350, 400, 446, 450, 452, 482]
    TRP_6_SIGNATURE = [19, 37, 45, 73, 75, 90, 129, 130, 142, 157, 181, 192, 194, 219,
                       221, 225, 227, 237, 287, 306, 337, 339, 350, 353, 356, 411, 462, 505]

    TRP_5_SIGNATURE_RESIDUES = "VSILIREPGLPRGVPRAVLPGEA"
    TRP_6_SIGNATURE_RESIDUES = "TEGCAGFDAYHDRFGNADYGLSIIAKIL"

    @staticmethod
    def get_signatures() -> List[List[int]]:
        raise NotImplementedError()


class TryptophanSubstrate(FlavinDependents):
    def __init__(self, name: str, residues: str) -> None:
        super().__init__(name, residues)
        self.cutoff = None

    def update_match(self, halogenase: HalogenaseResult, hit: HalogenaseHmmResult) -> None:
        if self.name == "trp_5":
            self.check_for_match(halogenase, hit, 5,
                                 cutoffs=[350, 850],
                                 sig_residues=FlavinDependents.TRP_5_SIGNATURE_RESIDUES)
            halogenase.substrate = "tryptophan"
        elif self.name == "trp_6_7":
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
        if hit.query_id == "trp_5":
            residues = get_residues(cds.translation, hit, [FlavinDependents.TRP_5_SIGNATURE])
        if hit.query_id == "trp_6_7":
            residues = get_residues(cds.translation, hit, [FlavinDependents.TRP_6_SIGNATURE])
        return residues


def run_halogenase_phmms(cluster_fasta: str) -> dict[str, list[HalogenaseHmmResult]]:
    """ Check if protein sequences hit any pHMM

        Arguments:
            cluster_fasta: string of protein sequences in a fasta format

        Returns:
            if there is a hit, it returns the properties of that hit
    """
    halogenase_hmms_by_id: dict = defaultdict(list)
    for sig in pHMM_SIGNATURES:
        sig.path = path.get_full_path(f'{sig.hmm_file}/{sig.name}')

        runresults = subprocessing.run_hmmsearch(sig.path, cluster_fasta)
        for runresult in runresults:
            for hsp in runresult.hsps:
                if hsp.bitscore > sig.cutoff:
                    hit = HalogenaseHmmResult(hsp.hit_id, hsp.bitscore,
                                              hsp.query_id, sig.name, sig.path)
                    halogenase_hmms_by_id[hsp.hit_id].append(hit)

    return halogenase_hmms_by_id


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

    results = subprocessing.hmmpfam.run_hmmpfam2(hmm_result.profile, f">query\n{sequence}", extra_args=args)
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

    sites = utils.extract_by_reference_positions(query, profile, [p - offset for p in positions if offset < p])

    return sites


def get_residues(translation: str, hmm_result: HalogenaseHmmResult,
                 signatures: list[list[int]]) -> dict[str, Optional[str]]:
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

    for signature in signatures:

        residue = search_signature_residues(translation, signature, hmm_result)
        signature_residues[hmm_result.query_id] = residue

    return signature_residues


def check_for_fdh(cds: CDSFeature, halogenase: HalogenaseResult, hit: HalogenaseHmmResult) -> None:
    """ Check if protein could be categorized as a Flavin-dependent enzyme

        Arguments:
            cds: gene/CDS and its properties
            halogenase: enzyme categorized as halogenase
            hit: details of the hit (e.g. bitscore, name of the profile, etc.)

        Returns:
            Modifies the HalogenasesResults by adding the label of the enzyme family,
            and possibly the signature and position of halogenation determined.
            If it doesn't meet the requirements, it doesn't make changes on the original instance.
    """
    signature_residues = TryptophanSubstrate.get_residues(cds, hit)
    if signature_residues:
        specific_signature_residues = signature_residues[hit.query_id]
        if not specific_signature_residues:
            return
        halogenase_category = TryptophanSubstrate(hit.query_id, specific_signature_residues)
        halogenase_category.update_match(halogenase, hit)


def check_for_halogenases(cds: CDSFeature, hmm_results: list[HalogenaseHmmResult]
                          ) -> Optional[HalogenaseResult]:
    """ Categorizes enzymes based on wether they hit the pHMM
        and have the required signature residues

        Arguments:
            cds: the query cds
            hmm_results: list of the halogenase hits and their properties
        Returns:
            categorized halogenase
    """
    if not hmm_results:
        logging.debug("Hmmsearch did not return any hit.")
        return None

    halogenase_match = HalogenaseResult(cds.get_name())
    for hit in hmm_results:
        check_for_fdh(cds, halogenase_match, hit)

    return halogenase_match


def specific_analysis(record: Record) -> list[HalogenaseResult]:
    """ Categorization of enzyme, categorizes any halogenase in a cds in regions

        Arguments: record instance,
                   which holds information of the identified clusters
        Returns: list of HalogenasesResults instances representing halogenase enzymes,
                 if there is a clear best match for a given enzyme, then the information
                 about the position of halogenation, the confidence of the categorization,
                 and the characteristic residues is provided.
                 If the enzyme can be categorized into several groups, with the same confidence,
                 then the above mentioned informations are not defined,
                 and the information about the catogries is in the potential_enzymes attribute.
    """
    potential_enzymes: list[HalogenaseResult] = []

    features = record.get_cds_features_within_regions()
    hmmsearch_fasta = fasta.get_fasta_from_features(features)

    hmm_hits = run_halogenase_phmms(hmmsearch_fasta)

    for protein, hits in hmm_hits.items():
        cds = record.get_cds_by_name(protein)
        potential_enzyme = check_for_halogenases(cds, hits)
        if potential_enzyme:
            potential_enzymes.append(potential_enzyme)

    for enzyme in potential_enzymes:
        enzyme.finalize_enzyme()

    return potential_enzymes
