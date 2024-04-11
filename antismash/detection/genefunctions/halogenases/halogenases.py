# License: GNU Affero General Public License v3 or later
# A copy of GNU AGPL v3 should have been included in this software package in LICENSE.txt.

# for test files, silence irrelevant and noisy pylint warnings
# pylint: disable=use-implicit-booleaness-not-comparison,protected-access,missing-docstring

from typing import Any, List, Optional, Iterable, Union, Dict

from dataclasses import dataclass, field

from antismash.common.hmmscan_refinement import HMMResult

NAME = "halogenases_analysis"

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
    cofactor: str
    family: str
    confidence: float
    signature: Union[str, Dict[str,str]]
    substrate: Optional[Union[int, str]] = None
    position: Optional[Union[int, str, List[int]]] = None
    number_of_decorations: Optional[dict[str, int]] = None

    def to_json(self) -> dict[str, Any]:
        return vars(self)

    @classmethod
    def from_json(cls, data: dict[str, Any]) -> "Match":
        return cls(**data)

@dataclass
class FlavinDependentHalogenases:
    cds_name: str
    cofactor: str
    family: str
    substrates: Union[str, List[str], None] = None
    target_positions: Optional[dict[str, int]] = None
    number_of_decorations: Optional[dict[str, int]] = None
    consensus_residues: Optional[dict[str, str]] = None
    confidence: float = 0
    potential_matches: list[Match] = field(default_factory=list)

    def add_potential_matches(self, match: Match) -> None:
        """ Adds the features of an enzyme group to list"""
        self.potential_matches.append(match)

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
                if abs(profile.confidence-highest_confidence) <= 0.005:
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
            self.cofactor = best_match.cofactor
            self.family = best_match.family
            self.target_positions = best_match.position
            self.consensus_residues = best_match.signature
            self.confidence = best_match.confidence

            if best_match.substrate:
                self.substrates = best_match.substrate

    def to_json(self) -> Dict[str, Any]:
        """ Constructs a JSON representation of this instance """
        potential_matches_json = [match.to_json() for match in self.potential_matches]

        return {
            "cds_name": self.cds_name,
            "family": self.family,
            "cofactor": self.cofactor,
            "substrates": self.substrates,
            "target_positions": self.target_positions,
            "number_of_decorations": self.number_of_decorations,
            "consensus_residues": self.consensus_residues,
            "confidence": self.confidence,
            "potential_matches": potential_matches_json
        }

    @classmethod
    def from_json(cls, data: Dict[str, Any]) -> "FlavinDependentHalogenases":
        """ Constructs the TailoringEnzymes from the JSON representation """

        cds_name = data["cds_name"]
        family = data["family"]
        cofactor = data["cofactor"]
        substrates = data["substrates"]
        target_positions = data["target_positions"]
        number_of_decorations = data["number_of_decorations"]
        consensus_residues = data["consensus_residues"]
        confidence = data["confidence"]
        potential_matches = [Match.from_json(profile) for profile in data["potential_matches"]]
        enzyme = cls(cds_name, cofactor, family, substrates, target_positions,
                     number_of_decorations, consensus_residues, confidence, potential_matches)
        return enzyme

class TailoringEnzymes():
    init: str
