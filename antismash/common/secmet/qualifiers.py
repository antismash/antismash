# License: GNU Affero General Public License v3 or later
# A copy of GNU AGPL v3 should have been included in this software package in LICENSE.txt.

""" Classes representing complex qualifiers for features.
"""

import bisect
from collections import defaultdict
from enum import Enum, unique
import re
from typing import Any, Dict, Iterable, Iterator, List, Set, Sequence, Tuple, Optional, Union  # pylint: disable=unused-import


def parse_format(fmt: str, data: str) -> Sequence[str]:
    """ Reverse of str.format(), pulls values from an input string that match
        positions of {} in a format string. Raises a ValueError if the match
        cannot be found.
    """
    safe = fmt.replace('(', r'\(').replace(')', r'\)')
    regex = "^{}$".format(safe.replace("{}", "(.+?)"))
    res = re.search(regex, data)
    if res is None:
        raise ValueError("Could not match format %r to input %r" % (fmt, data))
    return res.groups()


class ActiveSiteFinderQualifier:
    """ A qualifier for tracking active sites found (or not found) within a CDS.
    """
    def __init__(self) -> None:
        self._hits = set()  # type: Set[str]

    @property
    def hits(self) -> List[str]:
        """ Returns a list of all active site notes """
        return sorted(list(self._hits))

    def add(self, label: str) -> None:
        """ Adds an active site presence label to the qualifier. """
        self._hits.add(str(label))

    def to_biopython(self) -> List[str]:
        """ Creates a BioPython style qualifier from the qualifier """
        return self.hits

    def __bool__(self) -> bool:
        return bool(self._hits)


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
            if feature_name:
                self.feature_name = str(feature_name)
            else:
                self.feature_name = None
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
        self.predictions = {}  # type: Dict[str, str]
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

        def to_json(self) -> List[Union[str, float, int]]:
            """ Constructs a JSON-friendly representation of a Domain """
            return [self.query_id, self.evalue, self.bitscore, self.nseeds, self.tool]

        @classmethod
        def from_string(cls, line: str) -> "SecMetQualifier.Domain":
            """ Rebuilds a Domain from a string (e.g. from a genbank file) """
            return cls.from_json(parse_format(cls.qualifier_label, line))

        @classmethod
        def from_json(cls, json: Sequence[Union[str, float]]) -> "SecMetQualifier.Domain":
            """ Rebuilds a Domain from a JSON representation """
            assert len(json) == 5, json
            return cls(str(json[0]), float(json[1]), float(json[2]), int(json[3]), str(json[4]))

    def __init__(self, products: Set[str], domains: List["SecMetQualifier.Domain"]) -> None:
        self._domains = domains
        self.domain_ids = []  # type: List[str]
        for domain in self._domains:
            assert isinstance(domain, SecMetQualifier.Domain)
            self.domain_ids.append(domain.query_id)
        self._products = set()  # type: Set[str]
        self.add_products(products)
        self.kind = "biosynthetic"
        super().__init__()

    def __iter__(self) -> Iterator[str]:
        yield "Type: %s" % self.clustertype
        yield "; ".join(map(str, self._domains))
        yield "Kind: %s" % self.kind

    def append(self, _item: Any) -> None:
        raise NotImplementedError("Appending to this list won't work")

    def extend(self, _items: Iterable[Any]) -> None:
        raise NotImplementedError("Extending this list won't work")

    def add_products(self, products: Set[str]) -> None:
        """ Adds one or more products to the qualifier """
        assert isinstance(products, set), type(products)
        for product in products:
            assert isinstance(product, str) and "-" not in product, product
        self._products.update(products)

    def add_domains(self, domains: List["SecMetQualifier.Domain"]) -> None:
        """ Add a group of Domains to the the qualifier """
        for domain in domains:
            assert isinstance(domain, SecMetQualifier.Domain)
        self._domains.extend(domains)

    @property
    def domains(self) -> List["SecMetQualifier.Domain"]:
        """ A list of domains stored in the qualifier"""
        return list(self._domains)

    @property
    def products(self) -> List[str]:
        """ A list of all products a feature is involved in"""
        return sorted(self._products)

    @property
    def clustertype(self) -> str:
        """ A string hypen-separated products """
        return "-".join(sorted(self.products))

    @staticmethod
    def from_biopython(qualifier: List[str]) -> "SecMetQualifier":
        """ Converts a BioPython style qualifier into a SecMetQualifier. """
        domains = []
        products = set()  # type: Set[str]
        kind = "biosynthetic"
        if len(qualifier) != 3:
            raise ValueError("Cannot parse qualifier: %s" % qualifier)
        for value in qualifier:
            if value.startswith("Type: "):
                products = set(value.split("Type: ", 1)[1].split("-"))
            elif value.startswith("Kind: "):
                kind = value.split("Kind: ", 1)[1]
                assert kind == "biosynthetic", kind  # since it's the only kind we have
            else:
                domain_strings = value.split("; ")
                for domain_string in domain_strings:
                    domains.append(SecMetQualifier.Domain.from_string(domain_string))
        if not (domains and products and kind):
            raise ValueError("Cannot parse qualifier: %s" % qualifier)
        return SecMetQualifier(products, domains)

    def __len__(self) -> int:
        return 3


@unique
class GeneFunction(Enum):
    """ An Enum representing the function of a gene.
        Allows for more flexible conversion and more robust value constraints.
    """
    OTHER = 0
    CORE = 1
    ADDITIONAL = 2
    TRANSPORT = 3
    REGULATORY = 4

    def __str__(self) -> str:
        # because this information ends up in the record, make these more
        # labels more meaningful for users
        if self == GeneFunction.CORE:
            return "biosynthetic"
        if self == GeneFunction.ADDITIONAL:
            return "biosynthetic-additional"
        return str(self.name).lower()

    @staticmethod
    def from_string(label: str) -> "GeneFunction":
        """ Converts a string to a GeneFunction instance when possible.
            Raises an error if not possible.
        """
        for value in GeneFunction:
            if str(value) == label:
                return value
        raise ValueError("Unknown gene function label: %s" % label)


class GeneFunctionAnnotations:
    """ Tracks multiple annotations of gene functions. Each annotation contains
        the declared function of the gene, which tool determined the function,
        and a comment on how the determination was made.
    """
    slots = ["_annotations", "_by_tool", "_by_function"]

    class GeneFunctionAnnotation:  # pylint: disable=too-few-public-methods
        """ A single instance of an annotation. """
        slots = ["function", "tool", "description"]

        def __init__(self, function: GeneFunction, tool: str, description: str) -> None:
            assert isinstance(function, GeneFunction), "wrong type: %s" % type(function)
            assert tool and len(tool.split()) == 1, tool  # no whitespace allowed in tool name
            assert description
            self.function = function
            self.tool = str(tool)
            self.description = str(description)

        def __str__(self) -> str:
            return "%s (%s) %s" % (self.function, self.tool, self.description)

        def __repr__(self) -> str:
            return "GeneFunctionAnnotation(function=%r, tool='%s', '%s')" % (self.function, self.tool, self.description)

    def __init__(self) -> None:
        self._annotations = []  # type: List["GeneFunctionAnnotations.GeneFunctionAnnotation"]
        self._by_tool = defaultdict(list)  # type: Dict[str, List["GeneFunctionAnnotations.GeneFunctionAnnotation"]]
        self._by_function = defaultdict(list)  # type: Dict[GeneFunction, List["GeneFunctionAnnotations.GeneFunctionAnnotation"]]

    def __iter__(self) -> Iterator["GeneFunctionAnnotations.GeneFunctionAnnotation"]:
        for annotation in self._annotations:
            yield annotation

    def __len__(self) -> int:
        return len(self._annotations)

    def add(self, function: GeneFunction, tool: str, description: str
            ) -> "GeneFunctionAnnotations.GeneFunctionAnnotation":
        """ Adds a gene function annotation. If an existing annotation has all
            the same values as those provided, a duplicate will not be added.

            Arguments:
                function: a GeneFunction value
                tool: a string naming the tool that determined the gene function
                description: description text for storing extra information

            Returns:
                the GeneFunction added (or the existing one if duplicated)
        """
        tool = str(tool)
        description = str(description)
        # if there's already an exactly similar function annotation, skip adding
        existing_functions = self._by_function.get(function, [])
        for existing in existing_functions:
            if existing.tool == tool and existing.description == description:
                return existing
        new = GeneFunctionAnnotations.GeneFunctionAnnotation(function, tool, description)
        self._by_tool[tool].append(new)
        self._by_function[function].append(new)
        self._annotations.append(new)
        return new

    def add_from_qualifier(self, qualifier: List[str]) -> None:
        """ Converts a string-based qualifier into one or more GeneFunctions and
            adds them to the current set.
            Expected format will be as per str(GeneFunctionAnnotation).
        """
        for section in qualifier:
            function, tool, description = section.split(maxsplit=2)
            tool = tool[1:-1]  # strip the ()
            self.add(GeneFunction.from_string(function), tool, description)

    def get_by_tool(self, tool: str) -> Optional[List]:
        """ Returns a list of all GeneFunctionAnnotations which were added with
            the tool name provided.
        """
        return self._by_tool.get(tool)

    def get_by_function(self, function: GeneFunction) -> Optional[List]:
        """ Returns a list of GeneFunctionAnnotations which have the same
            gene function as that provided.
        """
        return self._by_function.get(function)

    def get_classification(self) -> GeneFunction:
        """ Returns the function of the gene, in the priority order of priority:
             - CORE if cluster defined by this gene
             - the function determined by smCOGs, if it exists
             - a function from all tools if they all agree
             - or OTHER
        """
        # if no annotations, skip to OTHER
        if not self._annotations:
            return GeneFunction.OTHER
        # if any CORE function set, use that
        if self._by_function.get(GeneFunction.CORE):
            return GeneFunction.CORE
        # then priority for smcogs
        annotations = self._by_tool.get("smcogs")
        if annotations:
            return annotations[0].function
        # otherwise check all agree
        function = self._annotations[0].function
        for annotation in self._annotations[1:]:
            if annotation.function != function:
                return GeneFunction.OTHER
        return function

    def clear(self) -> None:
        """ Removes all gene functions from the annotation """
        self._annotations = []
        self._by_tool = defaultdict(list)
        self._by_function = defaultdict(list)


class GOQualifier:
    """A qualifier for tracking Gene Ontology terms for a PFAM domain.
        Cannot be directly used as a qualifier for BioPython's SeqFeature.
    """
    def __init__(self, go_entries: Dict[str, str]) -> None:  # dict mapping Gene Ontology IDs to readable descriptions
        self.go_entries = go_entries
        self.ids = list(go_entries)
        self.descriptions = list(go_entries.values())

    def to_biopython(self) -> List[str]:
        """Convert GOQualifier to BioPython-style qualifier."""
        return ["{}: {}".format(go_id, go_description) for go_id, go_description in sorted(self.go_entries.items())]

    @staticmethod
    def from_biopython(qualifier: List[str]) -> "GOQualifier":
        """Convert BioPython-style qualifier to GOQualifier.

            Arguments:
                qualifier: BioPython-style qualifier (list of strings)

            Returns:
                A GOQualifier instance constructed from the qualifier.

        """
        go_entries = {}
        for go_string in qualifier:
            go_id, separator, go_description = go_string.partition(": ")
            if not separator:
                raise ValueError("Cannot parse qualifier: %s" % qualifier)
            go_entries[go_id] = go_description
        result = GOQualifier(go_entries)
        return result
