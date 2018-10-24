# License: GNU Affero General Public License v3 or later
# A copy of GNU AGPL v3 should have been included in this software package in LICENSE.txt.

""" Annotations for secondary metabolites """

import re
from typing import Any, Iterable, Iterator, List, Sequence, Set, Union


def _parse_format(fmt: str, data: str) -> Sequence[str]:
    """ Reverse of str.format(), pulls values from an input string that match
        positions of {} in a format string. Raises a ValueError if no matches
        are found.

        If the {} sections have a specific formatter, those specifics will be ignored

        Arguments:
            fmt: the format to parse with
            data: the string that was formatted with fmt and specific values

        Returns:
            a sequence of strings, one for each value found
    """
    # escape anything that might cause issues when using it as a regex later
    safe = re.escape(fmt)
    # the simple search here would be {.*?} for a non-greedy brace pair, but
    # because python format strings use {{ and }} as literal braces, they need to be excluded
    # e.g. "{{{{{:g}}}}}".format(1e-5) == "{{1e-5}}"
    # so, step 1: unescape all braces
    safe = safe.replace(r"\{", "{").replace(r"\}", "}")
    # step 2: combine all double braces to an escaped brace
    safe = safe.replace("{{", r"\{")
    # but from the outside in, so reverse the search and reverse the replacement
    safe = safe[::-1].replace("}}", "}\\")[::-1]  # not a raw string because that breaks the interpreter
    # step 3: replace all unescaped brace pairs with a capture group
    sub_search = r"(?<!\\)({.*?(?<!\\)})"
    # core    "({.*?})" any brace pair and its contents
    # special "(?<!\\)" disallows any core starting or ending with an escaped brace
    # this combo allows capturing the nested parts correctly
    regex = "^{}$".format(re.sub(sub_search, r"(.+?)", safe))

    res = re.search(regex, data)
    if res is None:
        raise ValueError("could not match format %r to input %r" % (fmt, data))
    # don't return matches, since that includes the braces, just use the groups
    return res.groups()


class SecMetQualifier(list):
    """ A qualifier for tracking various secondary metabolite information about
        a CDS.

        Can be directly used as a qualifier for BioPython's SeqFeature.
    """
    class Domain:
        """ A simple container for the information needed to create a domain """
        qualifier_label = "{} (E-value: {}, bitscore: {}, seeds: {}, tool: {})"

        def __init__(self, name: str, evalue: float, bitscore: float, nseeds: int,
                     tool: str) -> None:
            self.query_id = str(name)
            self.evalue = float(evalue)
            self.bitscore = float(bitscore)
            self.nseeds = int(nseeds)
            self.tool = str(tool)

        def __repr__(self) -> str:
            return str(self)

        def __str__(self) -> str:
            return self.qualifier_label.format(self.query_id, self.evalue,
                                               self.bitscore, self.nseeds, self.tool)

        def __eq__(self, other: Any) -> bool:
            if not isinstance(other, type(self)):
                return False
            return (self.query_id == other.query_id
                    and self.evalue == other.evalue
                    and self.bitscore == other.bitscore
                    and self.nseeds == other.nseeds
                    and self.tool == other.tool)

        def to_json(self) -> List[Union[str, float, int]]:
            """ Constructs a JSON-friendly representation of a Domain """
            return [self.query_id, self.evalue, self.bitscore, self.nseeds, self.tool]

        @classmethod
        def from_string(cls, line: str) -> "SecMetQualifier.Domain":
            """ Rebuilds a Domain from a string (e.g. from a genbank file) """
            return cls.from_json(_parse_format(cls.qualifier_label, line))

        @classmethod
        def from_json(cls, json: Sequence[Union[str, float]]) -> "SecMetQualifier.Domain":
            """ Rebuilds a Domain from a JSON representation """
            assert len(json) == 5, json
            return cls(str(json[0]), float(json[1]), float(json[2]), int(json[3]), str(json[4]))

    def __init__(self, domains: List["SecMetQualifier.Domain"]) -> None:
        self._domains = []  # type: List["SecMetQualifier.Domain"]
        self.domain_ids = []  # type: List[str]
        self.unique_domain_ids = set()  # type: Set[str]
        self.add_domains(domains)
        self.kind = "biosynthetic"
        super().__init__()

    def __iter__(self) -> Iterator[str]:
        yield "; ".join(map(str, self._domains))
        yield "Kind: %s" % self.kind

    def __getitem__(self, selection: Union[slice, int]) -> str:  # type: ignore
        return str(list(self)[selection])

    def append(self, _item: Any) -> None:
        raise NotImplementedError("Appending to this list won't work")

    def extend(self, _items: Iterable[Any]) -> None:
        raise NotImplementedError("Extending this list won't work")

    def add_domains(self, domains: List["SecMetQualifier.Domain"]) -> None:
        """ Add a group of Domains to the the qualifier """
        unique = []
        for domain in domains:
            assert isinstance(domain, SecMetQualifier.Domain)
            if domain.query_id in self.unique_domain_ids:
                continue  # no sense keeping duplicates
            self.unique_domain_ids.add(domain.query_id)
            unique.append(domain)
            self.domain_ids.append(domain.query_id)
        self._domains.extend(unique)

    @property
    def domains(self) -> List["SecMetQualifier.Domain"]:
        """ A list of domains stored in the qualifier"""
        return list(self._domains)

    @staticmethod
    def from_biopython(qualifier: List[str]) -> "SecMetQualifier":
        """ Converts a BioPython style qualifier into a SecMetQualifier. """
        domains = []
        kind = "biosynthetic"
        if len(qualifier) != 2:
            raise ValueError("Cannot parse qualifier: %s" % qualifier)
        for value in qualifier:
            if value.startswith("Kind: "):
                kind = value.split("Kind: ", 1)[1]
                assert kind == "biosynthetic", kind  # since it's the only kind we have
            else:
                domain_strings = value.split("; ")
                for domain_string in domain_strings:
                    domains.append(SecMetQualifier.Domain.from_string(domain_string))
        if not (domains and kind):
            raise ValueError("Cannot parse qualifier: %s" % qualifier)
        return SecMetQualifier(domains)

    def __len__(self) -> int:
        return 2
