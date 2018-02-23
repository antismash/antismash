# License: GNU Affero General Public License v3 or later
# A copy of GNU AGPL v3 should have been included in this software package in LICENSE.txt.

""" Classes representing complex qualifiers for features.
"""

from collections import defaultdict
from enum import Enum, unique
from typing import Dict, List, Set, Optional, Union


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
            self.label = str(label)
            self.name = str(name)
            self.start = int(start)
            self.end = int(end)
            self.evalue = float(evalue)
            self.bitscore = float(bitscore)
            self.predictions = {}  # type: Dict[str, str] # method to prediction name

        def __repr__(self):
            return "NRPSPKSQualifier.Domain(%s, label=%s, start=%d, end=%d)" % (
                        self.name, self.label, self.start, self.end)

    def __init__(self) -> None:
        super().__init__()
        self.type = "uninitialised"
        self.subtypes = []  # type: List[str]
        self.domains = []  # type: List["NRPSPKSQualifier.Domain"]
        self.domain_names = []  # type: List[str]
        self.predictions = {}
        self.cal_counter = 0
        self.at_counter = 0
        self.kr_counter = 0
        self.a_counter = 0
        self.ks_counter = 0
        self.other_counter = 0

    def append(self, _value):
        raise NotImplementedError("Appending to this list won't work, use add_subtype() or add_domain()")

    def extend(self, _values):
        raise NotImplementedError("Extending this list won't work")

    def __len__(self):
        return len(self.subtypes) + len(self.domains)

    def __iter__(self):
        for domain in self.domains:
            yield "NRPS/PKS Domain: %s (%s-%s). E-value: %s. Score: %s;" % (domain.name,
                    domain.start, domain.end, domain.evalue, domain.bitscore)
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

        self.domains.append(NRPSPKSQualifier.Domain(domain.hit_id, suffix,
                domain.query_start, domain.query_end, domain.evalue, domain.bitscore))
        self.domain_names.append(domain.hit_id)


class SecMetQualifier(list):
    """ A qualifier for tracking various secondary metabolite information about
        a CDS.

        Can be directly used as a qualifier for BioPython's SeqFeature.
    """
    class Domain:
        """ A simple container for the information needed to create a domain """
        def __init__(self, name: str, evalue: float, bitscore: float, nseeds: str,
                     tool: str) -> None:
            self.query_id = str(name)
            self.evalue = float(evalue)
            self.bitscore = float(bitscore)
            self.nseeds = str(nseeds)  # TODO: will be int once all HMM inputs are cleaned up
            self.tool = str(tool)

        def __repr__(self):
            return str(self)

        def __str__(self):
            ret = "{} E-value: {}, bitscore: {}, seeds: {}"
            return ret.format(self.query_id, self.evalue, self.bitscore, self.nseeds)

    def __init__(self, products: Set[str], domains: Union[List["SecMetQualifier.Domain"], List[str]]) -> None:
        self._domains = domains  # Domain instance or str
        self.domain_ids = []  # type: List[str]
        if domains and not isinstance(domains[0], str):  # SecMetResult
            for domain in self._domains:
                assert isinstance(domain, SecMetQualifier.Domain)
                self.domain_ids.append(domain.query_id)
        else:  # TODO: regenerate a Domain from the string
            for domain in self._domains:
                assert isinstance(domain, str)
                self.domain_ids.append(domain.split()[0])
        self._products = set()  # type: Set[str]
        self.add_products(products)
        self.kind = "biosynthetic"
        super().__init__()

    def __iter__(self):
        yield "Type: %s" % self.clustertype
        yield "; ".join(map(str, self._domains))
        yield "Kind: %s" % self.kind

    def append(self, _item):
        raise NotImplementedError("Appending to this list won't work")

    def extend(self, _items):
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
    def from_biopython(qualifier) -> "SecMetQualifier":
        """ Converts a BioPython style qualifier into a SecMetQualifier. """
        domains = None
        products = None
        kind = None
        if len(qualifier) != 3:
            raise ValueError("Cannot parse qualifier: %s" % qualifier)
        for value in qualifier:
            if value.startswith("Type: "):
                products = set(value.split("Type: ", 1)[1].split("-"))
            elif value.startswith("Kind: "):
                kind = value.split("Kind: ", 1)[1]
                assert kind == "biosynthetic", kind  # since it's the only kind we have
            else:
                domains = value.split("; ")
        if not (domains and products and kind):
            raise ValueError("Cannot parse qualifier: %s" % qualifier)
        return SecMetQualifier(products, domains)

    def __len__(self):
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
        raise ValueError("Unknown gene function label: %s", label)


class GeneFunctionAnnotations:
    """ Tracks multiple annotations of gene functions. Each annotation contains
        the declared function of the gene, which tool determined the function,
        and a comment on how the determination was made.
    """
    slots = ["_annotations", "_by_tool", "_by_function"]

    class GeneFunctionAnnotation:
        """ A single instance of an annotation. """
        slots = ["function", "tool", "description"]

        def __init__(self, function: GeneFunction, tool: str, description: str) -> None:
            assert isinstance(function, GeneFunction), "wrong type: %s" % type(function)
            assert tool and len(tool.split()) == 1, tool  # no whitespace allowed in tool name
            assert description
            self.function = function
            self.tool = str(tool)
            self.description = str(description)

        def __str__(self):
            return "%s (%s) %s" % (self.function, self.tool, self.description)

        def __repr__(self):
            return "GeneFunctionAnnotation(function=%r, tool='%s', '%s')" % (self.function, self.tool, self.description)

    def __init__(self):
        self._annotations = []
        self._by_tool = defaultdict(list)
        self._by_function = defaultdict(list)

    def __iter__(self):
        for annotation in self._annotations:
            yield annotation

    def __len__(self):
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

    def add_from_qualifier(self, qualifier: List[str]):
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
