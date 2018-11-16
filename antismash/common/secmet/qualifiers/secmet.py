# License: GNU Affero General Public License v3 or later
# A copy of GNU AGPL v3 should have been included in this software package in LICENSE.txt.

""" Annotations for secondary metabolites """

import re
from typing import Any, Iterator, List, Sequence, Union
from typing import Set  # comment hints  # pylint: disable=unused-import


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


class SecMetQualifier:
    """ A qualifier for tracking various secondary metabolite information about
        a CDS.
    """
    class Domain:
        """ A simple container for the information needed to create a domain """
        qualifier_label = "{} (E-value: {}, bitscore: {}, seeds: {}, tool: {})"

        def __init__(self, name: str, evalue: float, bitscore: float, nseeds: int,
                     tool: str) -> None:
            self.name = str(name)
            self.evalue = float(evalue)
            self.bitscore = float(bitscore)
            self.nseeds = int(nseeds)
            self.tool = str(tool)

        def __repr__(self) -> str:
            return str(self)

        def __str__(self) -> str:
            return self.qualifier_label.format(self.name, self.evalue,
                                               self.bitscore, self.nseeds, self.tool)

        def __eq__(self, other: Any) -> bool:
            if not isinstance(other, type(self)):
                return False
            return (self.name == other.name
                    and self.evalue == other.evalue
                    and self.bitscore == other.bitscore
                    and self.nseeds == other.nseeds
                    and self.tool == other.tool)

        def to_json(self) -> List[Union[str, float, int]]:
            """ Constructs a JSON-friendly representation of a Domain """
            return [self.name, self.evalue, self.bitscore, self.nseeds, self.tool]

        @classmethod
        def from_string(cls, line: str) -> "SecMetQualifier.Domain":
            """ Rebuilds a Domain from a string (e.g. from a genbank file) """
            assert isinstance(line, str), type(line)
            return cls.from_json(_parse_format(cls.qualifier_label, line))

        @classmethod
        def from_json(cls, json: Sequence[Union[str, float]]) -> "SecMetQualifier.Domain":
            """ Rebuilds a Domain from a JSON representation """
            assert len(json) == 5, json
            return cls(str(json[0]), float(json[1]), float(json[2]), int(json[3]), str(json[4]))

    def __init__(self, domains: List["SecMetQualifier.Domain"] = None) -> None:
        self._domains = []  # type: List["SecMetQualifier.Domain"]
        self.domain_ids = []  # type: List[str]
        self.unique_domain_ids = set()  # type: Set[str]
        if domains is not None:
            self.add_domains(domains)
        super().__init__()

    def __iter__(self) -> Iterator["SecMetQualifier.Domain"]:
        return iter(self._domains)

    def add_domains(self, domains: List["SecMetQualifier.Domain"]) -> None:
        """ Add a group of Domains to the the qualifier """
        unique = []
        for domain in domains:
            assert isinstance(domain, SecMetQualifier.Domain)
            if domain.name in self.unique_domain_ids:
                continue  # no sense keeping duplicates
            self.unique_domain_ids.add(domain.name)
            unique.append(domain)
            self.domain_ids.append(domain.name)
        self._domains.extend(unique)

    @property
    def domains(self) -> List["SecMetQualifier.Domain"]:
        """ A list of domains stored in the qualifier"""
        return list(self._domains)

    @staticmethod
    def from_biopython(qualifier: List[str]) -> "SecMetQualifier":
        """ Converts a BioPython style qualifier into a SecMetQualifier. """
        domains = []
        for value in qualifier:
            domains.append(SecMetQualifier.Domain.from_string(value))
        return SecMetQualifier(domains)

    def __len__(self) -> int:
        return len(self._domains)
