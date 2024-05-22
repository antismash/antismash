# License: GNU Affero General Public License v3 or later
# A copy of GNU AGPL v3 should have been included in this software package in LICENSE.txt.

""" Annotations for NRPS/PKS domains """

import bisect
from dataclasses import dataclass, field
from typing import Any, Dict, Iterator, List, Tuple

from .secmet import _parse_format

_DOMAIN_FORMAT = "Domain: {} ({:d}-{:d}). E-value: {}. Score: {}. Matches aSDomain: {}"
_TYPE_FORMAT = "type: {}"


class _HMMResultLike:
    """ A class that has compatible members with antismash.common.hmmscan_refinement.HMMResult.
        Used for reconstructing domains from qualifiers
    """
    def __init__(self, hit_id: str, query_start: int, query_end: int,
                 evalue: float, bitscore: float, subtypes: List[str]) -> None:
        self.hit_id = hit_id
        self.query_start = query_start
        self.query_end = query_end
        self.evalue = evalue
        self.bitscore = bitscore
        self.detailed_names = subtypes


class NRPSPKSQualifier:
    """ A qualifier for tracking information about NRPS/PKS domains within a CDS.

        Can be used directly as a qualifier for Biopython's SeqFeature.
    """
    @dataclass
    class Domain:
        """ Contains information about a NRPS/PKS domain, including predictions
            made by modules.

            feature_name is identical to that of the AntismashDomain that contains
            this same information
        """
        name: str
        label: str
        start: int
        end: int
        evalue: float
        bitscore: float
        feature_name: str
        subtypes: list[str] = field(default_factory=list)
        _predictions: dict[str, str] = field(default_factory=dict)

        def __post_init__(self) -> None:
            if not self.feature_name:
                raise ValueError("a Domain must belong to a feature")

        def __lt__(self, other: "NRPSPKSQualifier.Domain") -> bool:
            return (self.start, self.end) < (other.start, other.end)

        def __repr__(self) -> str:
            return f"NRPSPKSQualifier.Domain({self.full_type}, label={self.label}, start={self.start}, end={self.end})"

        def get_predictions(self) -> Dict[str, str]:
            """ Returns a dictionary mapping prediction method to prediction """
            return dict(self._predictions)

        def add_prediction(self, method: str, prediction: str) -> None:
            """ Sets a prediction for a method to the domain """
            self._predictions[method] = prediction

        @property
        def full_type(self) -> str:
            """ Returns the type of a domain, including subtype, if present """
            if not self.subtypes:
                return self.name
            # limit to a single subtype, for the moment
            return f"{self.name}({')('.join(self.subtypes[:1])})"

    def __init__(self, strand: int) -> None:
        super().__init__()
        if strand not in [1, -1]:
            raise ValueError(f"strand must be 1 or -1, not {strand}")
        self.strand = strand
        self.type = "uninitialised"
        self._domains: List["NRPSPKSQualifier.Domain"] = []
        self._domain_names: List[str] = []
        self.cal_counter = 0
        self.at_counter = 0
        self.kr_counter = 0
        self.a_counter = 0
        self.ks_counter = 0
        self.other_counter = 0

    @property
    def domains(self) -> Tuple["NRPSPKSQualifier.Domain", ...]:
        """ Returns a list of Domains added to the qualifier, ordered by position """
        return tuple(self._domains)

    @property
    def domain_names(self) -> List[str]:
        """ Returns a list of domain names in order first to last position on the strand """
        return self._domain_names

    def __len__(self) -> int:
        return len(self._domains)

    def __iter__(self) -> Iterator[str]:
        for domain in self.domains:
            yield _DOMAIN_FORMAT.format(domain.full_type, domain.start, domain.end,
                                        domain.evalue, domain.bitscore, domain.feature_name)
        if self.type != "uninitialised":
            yield _TYPE_FORMAT.format(self.type)

    # the domain type Any is only to avoid circular dependencies
    def add_domain(self, domain: Any, feature_name: str) -> Domain:
        """ Adds a domain to the current set.

            Arguments:
                domain: the domain to add, this should be a HMMResult-like object
            (see: antismash.common.hmmscan_refinement.HMMResult).
                feature_name: the name of the matching AntismashDomain feature
                              in the same record as this qualifier

            Returns:
                None
        """
        assert not isinstance(domain, str)
        if domain.hit_id == "PKS_AT":
            self.at_counter += 1
            suffix = f"_AT{self.at_counter}"
        elif domain.hit_id == "PKS_KR":
            self.kr_counter += 1
            suffix = f"_KR{self.kr_counter}"
        elif domain.hit_id == "CAL_domain":
            self.cal_counter += 1
            suffix = f"_CAL{self.cal_counter}"
        elif domain.hit_id in ["AMP-binding", "A-OX"]:
            self.a_counter += 1
            suffix = f"_A{self.a_counter}"
        elif domain.hit_id == "PKS_KS":
            self.ks_counter += 1
            suffix = f"_KS{self.ks_counter}"
        else:
            self.other_counter += 1
            suffix = f"_OTHER{self.other_counter}"

        new = NRPSPKSQualifier.Domain(domain.hit_id, suffix,
                                      domain.query_start, domain.query_end,
                                      domain.evalue, domain.bitscore, feature_name,
                                      domain.detailed_names[1:])
        bisect.insort_right(self._domains, new)
        # update the domain name list
        self._domain_names = [domain.name for domain in self._domains]
        return new

    def add_from_qualifier(self, qualifiers: List[str]) -> None:
        """ Adds domains and types from a biopython-style qualifier list """
        if not qualifiers:
            return
        for qualifier in qualifiers:
            if qualifier.startswith("Domain: "):
                parts = _parse_format(_DOMAIN_FORMAT, qualifier)
                name = parts[0]
                types = []
                if "(" in name:
                    types = [chunk.strip(")") for chunk in parts[0].split("(")]
                    name = types[0]
                domain = _HMMResultLike(name, int(parts[1]), int(parts[2]),
                                        float(parts[3]), float(parts[4]), types)
                self.add_domain(domain, parts[5])
            elif qualifier.startswith("type: "):
                self.type = _parse_format(_TYPE_FORMAT, qualifier)[0]
            else:
                raise ValueError(f"unknown NRPS/PKS qualifier: {qualifier}")
