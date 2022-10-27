# License: GNU Affero General Public License v3 or later
# A copy of GNU AGPL v3 should have been included in this software package in LICENSE.txt.

""" Contains shared classes for the module.
"""

from dataclasses import dataclass
from typing import Any, Dict, List, Type, TypeVar

T = TypeVar("T", bound="JsonConvertible")
# groups of conserved amino acids for use in generating consensus strings
AMINO_GROUPS = {
    "A": "S",
    "D": "NE",
    "E": "QDK",
    "F": "YW",
    "H": "YN",
    "I": "VLM",
    "K": "RQE",
    "L": "MIV",
    "M": "LIV",
    "N": "DHS",
    "Q": "ERK",
    "R": "KQ",
    "S": "ANT",
    "T": "S",
    "V": "IML",
    "W": "YF",
    "Y": "FHW",
}


class JsonConvertible:
    """ A base class/mixin for all classes able to convert to and from JSON """
    def __init__(self, *, kwargs: Dict[str, Any]) -> None:
        raise NotImplementedError()

    def to_json(self) -> Dict[str, Any]:
        """ Converts the object into a JSON-compatible dictionary """
        return dict(vars(self))

    @classmethod
    def from_json(cls: Type[T], data: Dict[str, Any]) -> T:
        """ Rebuilds an instance from a JSON-compatible dictionary """
        return cls(**data)


@dataclass
class Segment(JsonConvertible):
    """ A class representing a sequence segment as reported by blast. The
        sequence covers from start position to end position, not the full length.
    """
    sequence: str
    start: int
    end: int
    full_length: int


@dataclass(frozen=True)
class Hit(JsonConvertible):
    """ A class containing all the relevant information for a hit.
    """
    BLAST_FORMAT = "6 qacc sacc nident qseq qstart qend qlen sseq sstart send slen"
    REFERENCE_ID_INDEX = 1
    query_name: str
    reference_fields: Dict[str, str]
    match_count: int
    query: Segment
    reference: Segment

    @property
    def similarity(self) -> float:
        """ Returns the similarity (between 0.0 and 1.0, inclusive) of the hit
            to the longest of query and reference sequences.
        """
        return self.match_count / max(self.query.full_length, self.reference.full_length)

    def get_consensus(self, space: str = " ") -> str:
        """ Returns the consensus string for the hit, with the same character for
            matches, '+' for strongly conserved mismatches, and the given "space"
            character for remaining mismatches.
        """
        return build_consensus_string(self.query.sequence, self.reference.sequence, space=space)

    @classmethod
    def from_simple_blast_line(cls, parts: List[str], fields: Dict[str, str]) -> "Hit":
        """ Constructs a Hit from a blast output line matching Hit.BLAST_FORMAT """
        return cls(
            parts[0],  # qacc
            fields,  # from sacc
            int(parts[2]),  # nident
            Segment(
                parts[3],  # qseq
                int(parts[4]),  # qstart
                int(parts[5]),  # qend
                int(parts[6]),  # qlen
            ),
            Segment(
                parts[7],  # sseq
                int(parts[8]),  # sstart
                int(parts[9]),  # send
                int(parts[10]),  # slen
            ),
        )

    def to_json(self) -> Dict[str, Any]:
        return {
            "query_name": self.query_name,
            "reference_fields": self.reference_fields,
            "query": self.query.to_json(),
            "reference": self.reference.to_json(),
            "match_count": self.match_count,
        }

    @classmethod
    def from_json(cls, data: Dict[str, Any]) -> "Hit":
        """ Rebuilds an instance from a JSON-compatible dictionary """
        for key in ("query", "reference"):
            data[key] = Segment.from_json(data[key])
        return cls(**data)


def build_consensus_string(first: str, second: str, space: str = " ") -> str:
    """ Returns the consensus string for the hit, with the same character for
        matches, '+' for strongly conserved mismatches, and the given "space"
        character for remaining mismatches.
    """
    assert len(first) == len(second)
    chars = []
    for query, reference in zip(first, second):
        if query == reference:
            chars.append(query)
            continue
        char = space
        if reference in AMINO_GROUPS.get(query, []):
            char = "+"
        chars.append(char)
    return "".join(chars)
