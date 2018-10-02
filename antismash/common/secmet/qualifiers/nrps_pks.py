# License: GNU Affero General Public License v3 or later
# A copy of GNU AGPL v3 should have been included in this software package in LICENSE.txt.

""" Annotations for NRPS/PKS domains """

import bisect
from typing import Any, Iterator, List, Tuple
from typing import Dict  # used in comment hints # pylint: disable=unused-import


class NRPSPKSQualifier(list):
    """ A qualifier for tracking information about NRPS/PKS domains within a CDS.

        Can be used directly as a qualifier for Biopython's SeqFeature.
    """
    class Domain:  # pylint: disable=too-few-public-methods
        """ Contains information about a NRPS/PKS domain, including predictions
            made by modules.

            feature_name is identical to that of the AntismashDomain that contains
            this same information
        """
        __slots__ = ["name", "label", "start", "end", "evalue", "bitscore",
                     "predictions", "feature_name"]

        def __init__(self, name: str, label: str, start: int, end: int,
                     evalue: float, bitscore: float, feature_name: str) -> None:
            self.label = str(label)
            self.name = str(name)
            self.start = int(start)
            self.end = int(end)
            self.evalue = float(evalue)
            self.bitscore = float(bitscore)
            if not feature_name:
                raise ValueError("a Domain must belong to a feature, feature_name is required")
            self.feature_name = str(feature_name)
            self.predictions = {}  # type: Dict[str, str] # method to prediction name

        def __lt__(self, other: "NRPSPKSQualifier.Domain") -> bool:
            return (self.start, self.end) < (other.start, other.end)

        def __repr__(self) -> str:
            return "NRPSPKSQualifier.Domain(%s, label=%s, start=%d, end=%d)" % (
                        self.name, self.label, self.start, self.end)

    def __init__(self, strand: int) -> None:
        super().__init__()
        if strand not in [1, -1]:
            raise ValueError("strand must be 1 or -1, not %s" % strand)
        self.strand = strand
        self.type = "uninitialised"
        self.subtypes = []  # type: List[str]
        self._domains = []  # type: List["NRPSPKSQualifier.Domain"]
        self._domain_names = []  # type: List[str]
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

    def append(self, _value: Any) -> None:
        raise NotImplementedError("Appending to this list won't work, use add_subtype() or add_domain()")

    def extend(self, _values: Any) -> None:
        raise NotImplementedError("Extending this list won't work")

    def __len__(self) -> int:
        return len(self.subtypes) + len(self._domains)

    def __iter__(self) -> Iterator[str]:
        for domain in self.domains:
            base = "NRPS/PKS Domain: %s (%s-%s). E-value: %s. Score: %s;"
            yield base % (domain.name, domain.start, domain.end, domain.evalue, domain.bitscore)
        for subtype in self.subtypes:
            yield "NRPS/PKS subtype: %s" % subtype

    def add_subtype(self, subtype: str) -> None:
        """ Adds a subtype to the existing list, e.g. 'Glycopeptide NRPS' or
            'NRPS-like protein'
        """
        assert isinstance(subtype, str)
        self.subtypes.append(subtype)

    # the domain type Any is only to avoid circular dependencies
    def add_domain(self, domain: Any, feature_name: str) -> None:
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
            suffix = "_AT%d" % self.at_counter
        elif domain.hit_id == "PKS_KR":
            self.kr_counter += 1
            suffix = "_KR%d" % self.kr_counter
        elif domain.hit_id == "CAL_domain":
            self.cal_counter += 1
            suffix = "_CAL%d" % self.cal_counter
        elif domain.hit_id == "AMP-binding":
            self.a_counter += 1
            suffix = "_A%d" % self.a_counter
        elif domain.hit_id == "PKS_KS":
            self.ks_counter += 1
            suffix = "_KS%d" % self.ks_counter
        else:
            self.other_counter += 1
            suffix = "_OTHER%d" % self.other_counter

        new = NRPSPKSQualifier.Domain(domain.hit_id, suffix,
                                      domain.query_start, domain.query_end,
                                      domain.evalue, domain.bitscore, feature_name)
        bisect.insort_right(self._domains, new)
        # update the domain name list
        self._domain_names = [domain.name for domain in self._domains]
