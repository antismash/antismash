# License: GNU Affero General Public License v3 or later
# A copy of GNU AGPL v3 should have been included in this software package in LICENSE.txt.

""" Data structures and related helpers for specific halogenase detection
    and classification
"""

from dataclasses import asdict, dataclass, field
from enum import StrEnum, auto
import os
from typing import Any, ClassVar, Iterator, Optional, Self

from antismash.common import path
from antismash.common.secmet import ECGroup, GeneFunction
from antismash.common.signature import HmmSignature

from ..core import FunctionResults, HMMHit

DEFAULT_PREFERENCE = 100

_BASE_DATA_PATH = path.get_full_path(__file__, "data")


class Conventionality(StrEnum):
    """ The options for a halogenases conventionality """
    AMBIGUOUS = auto()
    CONVENTIONAL = auto()
    UNCONVENTIONAL = auto()


CONVENTIONALITY_CONFIDENCES = {
    Conventionality.CONVENTIONAL: 1.,
    Conventionality.UNCONVENTIONAL: 0.,
}


def get_file_path(filename: str) -> str:
    """ Returns the path to a particular file, based on the location of this file

        Arguments:
            filename: the filename of the file

        Returns:
            a string with the combined path
    """
    return os.path.join(_BASE_DATA_PATH, filename)


@dataclass(slots=True)
class Match:  # pylint: disable=too-many-instance-attributes
    """ Contains details about which pHMM (profile) was hit, what position
        the halogenation occurs, what is the confidence of the categorization,
        and what are the signature residues of the protein sequence"""
    profile: str
    cofactor: str
    family: str
    confidence: float
    consensus_residues: str
    substrate: Optional[str] = None
    target_positions: Optional[tuple[int, ...]] = None
    number_of_decorations: str = ""

    @classmethod
    def from_json(cls, data: dict[str, Any]) -> "Match":
        """ Reconstructs an instance from JSON """
        data["target_positions"] = tuple(data["target_positions"])
        return cls(**data)


@dataclass(kw_only=True, slots=True)
class HalogenaseHit(HMMHit):
    """ A halogenase-specific HMM hit, allowing for internal HMM hits
    """
    profile: Optional["Profile"] = None
    conventionality: Conventionality = Conventionality.AMBIGUOUS
    conventionality_residues: dict[str, str] = field(default_factory=dict)
    potential_matches: list[Match] = field(default_factory=list)
    subfunctions: list[str] = field(default_factory=lambda: ["Halogenation"])

    def __post_init__(self) -> None:
        assert self.profile is not None

    @property
    def cofactor(self) -> str:
        """ The cofactor for the halogenase """
        assert self.profile
        return self.profile.cofactor

    @property
    def confidence(self) -> float:
        """ The confidence score of the hit """
        confidence = CONVENTIONALITY_CONFIDENCES.get(self.conventionality, 0.)

        if not self.potential_matches:
            return confidence

        # at least one match is good enough to bump up the confidence
        confidence += 1

        confidence += max(match.confidence for match in self.potential_matches)

        return confidence

    @property
    def family(self) -> str:
        """ The family of the relevant halogenase """
        assert self.profile
        return self.profile.family

    def get_best_matches(self, epsilon: float = 0.005) -> list[Match]:
        """ Gets the best motif matches for the instance.

            Arguments:
                epsilon: the allowable difference between the singular best hit and later hits

            Returns:
                a list of matches
        """
        best_matches: list[Match] = []
        if not self.potential_matches:
            return best_matches

        matches = sorted(self.potential_matches, key=lambda x: x.confidence, reverse=True)
        best_matches.append(matches[0])

        for match in best_matches[1:]:
            if best_matches[0].confidence - match.confidence <= epsilon:
                best_matches.append(match)

        return best_matches

    def get_full_description(self) -> str:
        if self.potential_matches:
            return f"{self.family} {self.reference_id.replace('_', ' ').replace('FDH', '').strip()} halogenase"
        return f"{self.conventionality} {self.family} halogenase"

    def is_conventional(self) -> bool:
        """ Returns True if the hit is explicitly a conventional halogenase """
        return self.conventionality == Conventionality.CONVENTIONAL

    @classmethod
    def from_json(cls, data: dict[str, Any]) -> Self:
        """ Constructs an instance from a JSON representation """
        potential_matches = [Match.from_json(profile) for profile in data.pop("potential_matches")]
        profile = Profile.from_json(data.pop("profile"))
        conventionality = Conventionality[data.pop("conventionality").upper()]
        return cls(**data, conventionality=conventionality, profile=profile, potential_matches=potential_matches)

    def to_json(self) -> dict[str, Any]:
        """ Constructs a JSON representation of this instance """
        data = asdict(self)
        assert self.profile
        data["profile"] = self.profile.to_json()
        return data


@dataclass(frozen=True, kw_only=True)
class MotifDetails:  # pylint: disable=too-many-instance-attributes
    """ A class for holding details about a specific motif.

        Attributes:
            name: the name of the motif
            positions: the positions of the motif's residues within a reference sequence
            residues: the residues at each of the provided positions
            substrate: the substrate to which this motif applies, if any
            decorations: a description of the decoration type or count for the motif
            maximum_distance: the maximum substitutions allowed for residue matching
            minimum_distance: the minimum substitutions required for residue matching
            target_positions: the positions at which the target will be modified
    """
    name: str
    positions: tuple[int, ...] = field(repr=False)
    residues: str
    substrate: str = ""
    decorations: str = ""
    preference_order: int = DEFAULT_PREFERENCE
    maximum_distance: int = 0
    minimum_distance: int = 0
    target_positions: tuple[int, ...] = field(default_factory=tuple)

    def __post_init__(self) -> None:
        assert not self.positions or isinstance(self.positions[0], int)
        assert len(list(self.positions)) == len(self.residues)
        assert tuple(sorted(self.positions)) == self.positions

    def __eq__(self, other: Any) -> bool:
        if isinstance(other, MotifDetails):
            return other.residues == self.residues
        if isinstance(other, str):
            return other == self.residues
        return False

    def __iter__(self) -> Iterator[tuple[int, str]]:
        return zip(self.positions, self.residues)

    def __len__(self) -> int:
        return len(self.residues)

    @classmethod
    def from_json(cls, data: dict[str, Any]) -> Self:
        """ Reconstruct an instance from JSON """
        data["positions"] = tuple(data["positions"])
        data["target_positions"] = tuple(data.get("target_positions", tuple()))
        return cls(**data)


@dataclass(frozen=True, kw_only=True)
class Profile:  # pylint: disable=too-many-instance-attributes
    """ A class for holding details about a particular halogenase type.

        Attributes:
            name: the identifier for the class
            description: a description of the halogenase type
            profile_name: the name of the matching pHMM
            filename: the name of the pHMM file itself
            cutoffs: a list of bitscore cutoffs ordered by increasing confidence
            cofactor: the cofactor of the particular halogenase
            family: the family of the particular halogenase
            motifs: details for each motif related to the halogenase type
            cutoff_step_multiplier: a penalty to that applies to any result using
                                    lower cutoffs, applied multiplicatively to each later cutoff
    """
    DEFAULT_PREFERENCE: ClassVar[int] = 100  # arbitrarily large position

    name: str
    description: str
    profile_name: str
    filename: str = field(repr=False)
    cutoffs: tuple[int, ...]

    cofactor: str = ""
    family: str = ""

    motifs: tuple[MotifDetails, ...] = field(default_factory=tuple)

    cutoff_step_multiplier: float = 0.5
    preference_penalty: float = 0.
    # some pre-cached values, since functools.cached_property ruins some documentation
    __motif_mapping: dict[str, MotifDetails] = field(repr=False, default_factory=dict)
    __hmm_profile: list[HmmSignature] = field(repr=False, default_factory=list)  # a list for mutability in init

    def __post_init__(self) -> None:
        # some workarounds for pre-caching derived properties within a frozen dataclass
        if sorted(self.cutoffs, reverse=True) != self.cutoffs:
            object.__setattr__(self, "cutoffs", sorted(self.cutoffs, reverse=True))

        self.__hmm_profile.append(HmmSignature(self.profile_name, self.description, self.cutoffs[-1], self.filename))
        self.__motif_mapping.update({motif.name: motif for motif in self.motifs})

        # value checks
        if len(self.__motif_mapping) != len(self.motifs):
            raise ValueError("provided motifs are not uniquely named")
        if not self.cutoffs:
            raise ValueError("at least one cutoff is required for the HMM profile")

        if not os.path.isabs(self.filename):
            object.__setattr__(self, "filename", get_file_path(self.filename))
            assert os.path.exists(self.filename)

    @property
    def profile(self) -> HmmSignature:
        """ The HMM profile for this halogenase type """
        return self.__hmm_profile[0]

    def get_matches_from_hit(self, retrieved_residues: dict[str, str], hit: HalogenaseHit,
                             confidence: float = 1.,
                             ) -> list[Match]:
        """ Searches for good motif matches for the halogenase type.

            Arguments:
                retrieved_residues: a dictionary mapping motif name to residues extracted for that motif
                hit: the HMMer hit from the halogenase type's profile
                confidence: a default confidence to use for all matches found

            Returns:
                a list of matches found
        """
        # if the hit is worse than the lowest of possible cutoffs, it's not relevant at all
        if hit.bitscore < self.cutoffs[-1]:
            return []

        # penalise the confidence if the highest cutoff not reached
        for cutoff in self.cutoffs:
            if hit.bitscore >= cutoff:
                break
            confidence *= self.cutoff_step_multiplier

        # no further processing required if the profile has not specific motifs to find
        if not self.motifs:
            return [self.create_match(confidence=confidence)]

        # gather all present motifs, ordered by preference order
        motifs_present = [(self.__motif_mapping[name], residues) for name, residues in retrieved_residues.items()]
        refined = []
        for motif, residues in motifs_present:
            # some motifs allow or require variation from the signature
            # if given in a regex-like form with `.` being ambiguous, that should count as equal
            distance = len(motif) - sum(q in (r, ".") for q, r in zip(residues, motif.residues))
            if motif.minimum_distance <= distance <= motif.maximum_distance:
                refined.append((motif, residues))
        motifs_present = sorted(refined, key=lambda motif_pair: motif_pair[0].preference_order)

        if not motifs_present:
            return []

        matches = []
        best_preference = motifs_present[0][0].preference_order
        for motif, residues in motifs_present:
            later_preference = motif.preference_order > best_preference
            if later_preference:
                confidence = max(0., confidence - self.preference_penalty)
            matches.append(self.create_match(confidence, residues, motif.substrate,
                                             motif.decorations, motif.target_positions))
        hit.potential_matches.extend(matches)
        return matches

    def create_match(self, confidence: float, residues: str = "", substrate: str = "",
                     decorations: str = "", target_positions: tuple[int, ...] = None,
                     ) -> Match:
        """ Creates a match instance for the halogenase type.

            Arguments:
                confidence: the confidence of the match
                residues: the residues of the match, as extracted for the motif
                substrate: the specific substrate of the match, if any
                decorations: decoration details, if relevant
                target_positions: the possible positions of substrate modification

            Returns:
                the new match
        """
        return Match(self.profile_name, self.cofactor, self.family,
                     confidence, residues,
                     target_positions=target_positions, substrate=substrate,
                     number_of_decorations=decorations)

    @property
    def motif_names(self) -> tuple[str, ...]:
        """ The names of the motifs within this halogenase type """
        return tuple(motif.name for motif in self.motifs)

    @classmethod
    def from_json(cls, data: dict[str, Any]) -> Self:
        """ Reconstructs an instance from JSON """
        data["motifs"] = [MotifDetails.from_json(motif) for motif in data["motifs"]]
        data["filename"] = get_file_path(data["filename"])
        return cls(**data)

    def to_json(self) -> dict[str, Any]:
        """ Converts the instance into JSON """
        return {key: val for key, val in asdict(self).items() if "__" not in key}


class HalogenaseResults(FunctionResults[HalogenaseHit]):
    """ A container for the results of halogenase classification """
    def __init__(self, best_hits: dict[str, HalogenaseHit] = None,
                 function_mapping: dict[str, GeneFunction] = None,
                 all_hits: dict[str, list[HalogenaseHit]] = None,
                 ) -> None:
        best_hits = best_hits or {}
        ec_groups = {name: [ECGroup.TRANSFERASES] for name in best_hits}
        subfunctions = {name: ["Halogenation"] for name, hit in best_hits.items()}
        super().__init__(tool="halogenases", best_hits=best_hits,
                         function_mapping=function_mapping or {},
                         group_mapping=ec_groups,
                         subfunction_mapping=subfunctions,
                         )
        self.all_hits = all_hits or {}

    @classmethod
    def from_json(cls, data: dict[str, Any]) -> Optional[Self]:
        mapping = {cds_name: GeneFunction.from_string(function)
                   for cds_name, function in data.pop("function_mapping").items()}
        hits = cls.regenerate_hits(data.pop("best_hits", {}))
        all_hits = {}
        for name, values in data.pop("all_hits", {}):
            all_hits[name] = [HalogenaseHit.from_json(value) for value in values]
        result = cls(best_hits=hits, function_mapping=mapping, all_hits=all_hits)
        return result

    @staticmethod
    def regenerate_hits(hits_by_name: dict[str, dict[str, Any]]) -> dict[str, HalogenaseHit]:
        return {name: HalogenaseHit.from_json(hit) for name, hit in hits_by_name.items()}


@dataclass
class Group:
    """ Contains a collection of profiles specific to a particular group of halogenases """
    profiles: dict[str, Profile]

    def get_matching_profiles(self, hit: HalogenaseHit) -> list[Profile]:
        """ Returns the matching profiles for a particular hit """
        return [profile for name, profile in self.profiles.items() if hit.reference_id == profile.profile_name]

    @classmethod
    def from_json(cls, data: dict[str, Any]) -> Self:
        """ Reconstructs an instance from JSON """
        profiles = [Profile.from_json(raw) for raw in data["profiles"]]
        return cls({profile.name: profile for profile in profiles})
