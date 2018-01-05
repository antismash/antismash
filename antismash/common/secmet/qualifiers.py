# License: GNU Affero General Public License v3 or later
# A copy of GNU AGPL v3 should have been included in this software package in LICENSE.txt.

""" Classes representing complex qualifiers for features.
"""

from typing import List


class NRPSPKSQualifier(list):
    """ A qualifier for tracking information about NRPS/PKS domains within a CDS.

        Can be used directly as a qualifier for Biopython's SeqFeature.
    """
    class Domain:
        """ Contains information about a NRPS/PKS domain, including predictions
            made by modules.
        """
        __slots__ = ["name", "label", "start", "end", "evalue", "bitscore",
                     "predictions"]

        def __init__(self, name: str, label: str, start: int, end: int,
                     evalue: float, bitscore: float) -> None:
            self.label = label
            self.name = name
            self.start = start
            self.end = end
            self.evalue = evalue
            self.bitscore = bitscore
            self.predictions = {}

    def __init__(self) -> None:
        super().__init__()
        self.type = "uninitialised"
        self.subtypes = []  # type: List[str]
        self.domains = []
        self.domain_names = []  # type: List[str]
        self.predictions = {}
        self.cal = 0
        self.at = 0
        self.kr = 0
        self.a = 0
        self.other = 0

    def append(self):
        raise NotImplementedError("Appending to this list won't work, use add_subtype() or add_domain()")

    def extend(self):
        raise NotImplementedError("Extending this list won't work")

    def __len__(self):
        return len(self.subtypes) + len(self.domains)

    def __iter__(self):
        for domain in self.domains:
            yield "NRPS/PKS Domain: %s (%s-%s). E-value: %s. Score: %s;" % (domain.hit_id,
                    domain.query_start, domain.query_end, domain.evalue, domain.bitscore)
        for subtype in self.subtypes:
            yield "NRPS/PKS subtype: %s" % subtype

    def add_subtype(self, subtype: str) -> None:
        """ Adds a subtype to the existing list, e.g. 'Glycopeptide NRPS' or
            'NRPS-like protein'
        """
        assert isinstance(subtype, str)
        self.subtypes.append(subtype)

    def add_domain(self, domain) -> None:
        """ Adds a domain to the current set.

            The domain should be a HMMResult-like object
            (see: antismash.common.hmmscan_refinement.HMMResult).
        """
        assert not isinstance(domain, str)
        if domain.hit_id == "PKS_AT":
            self.at += 1
            suffix = "_AT%d" % self.at
        elif domain.hit_id == "PKS_KR":
            self.kr += 1
            suffix = "_KR%d" % self.kr
        elif domain.hit_id == "CAL_domain":
            self.cal += 1
            suffix = "_CAL%d" % self.cal
        elif domain.hit_id == "AMP-binding":
            self.a += 1
            suffix = "_A%d" % self.a
        else:
            self.other += 1
            suffix = "_OTHER%d" % self.other

        self.domains.append(NRPSPKSQualifier.Domain(domain.hit_id, suffix,
                domain.query_start, domain.query_end, domain.evalue, domain.bitscore))
        self.domain_names.append(domain.hit_id)


class SecMetQualifier(list):
    """ A qualifier for tracking various secondary metabolite information about
        a CDS.

        Can be directly used as a qualifier for BioPython's SeqFeature.
    """
    def __init__(self, clustertype: str, domains: List) -> None:
        self.domains = domains  # SecMetResult instance or str
        if domains and not isinstance(domains[0], str):  # SecMetResult
            self.domain_ids = [domain.query_id for domain in self.domains]
        else:  # str
            self.domain_ids = [domain.split()[0] for domain in self.domains]
        self.clustertype = clustertype
        self.kind = "biosynthetic"
        super().__init__()

    def __iter__(self):
        yield "Type: %s" % self.clustertype
        yield "; ".join(map(str, self.domains))
        yield "Kind: %s" % self.kind

    def append(self, _item):
        raise NotImplementedError("Appending to this list won't work")

    def extend(self, _items):
        raise NotImplementedError("Extending this list won't work")

    @staticmethod
    def from_biopython(qualifier) -> "SecMetQualifier":
        """ Converts a BioPython style qualifier into a SecMetQualifier. """
        domains = None
        clustertype = None
        kind = None
        if len(qualifier) != 3:
            raise ValueError("Cannot parse qualifier: %s" % qualifier)
        for value in qualifier:
            if value.startswith("Type: "):
                clustertype = value.split("Type: ", 1)[0]
            elif value.startswith("Kind: "):
                kind = value.split("Kind: ", 1)[1]
                assert kind == "biosynthetic", kind  # since it's the only kind we have
            else:
                domains = value.split("; ")
        if not domains and clustertype and kind:
            raise ValueError("Cannot parse qualifier: %s" % qualifier)
        return SecMetQualifier(clustertype, domains)

    def __len__(self):
        return 3
