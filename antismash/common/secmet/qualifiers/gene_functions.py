# License: GNU Affero General Public License v3 or later
# A copy of GNU AGPL v3 should have been included in this software package in LICENSE.txt.

""" Annotations and helpers for annotating the function of genes """

from collections import defaultdict
from enum import Enum, unique
from typing import Any, Iterator, List, Optional
from typing import Dict  # comment hints  # pylint: disable=unused-import

from .secmet import _parse_format


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
    RESISTANCE = 5

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


class _GeneFunctionAnnotation:  # pylint: disable=too-few-public-methods
    """ A single instance of an annotation. """
    slots = ["function", "tool", "description", "product"]

    def __init__(self, function: GeneFunction, tool: str, description: str,
                 product: Optional[str]) -> None:
        assert isinstance(function, GeneFunction), "wrong type: %s" % type(function)
        assert tool and len(tool.split()) == 1, tool  # no whitespace allowed in tool name
        assert description
        if function == GeneFunction.CORE and not product:
            raise ValueError("CORE functions require a product to be specified")
        self.function = function
        self.tool = str(tool)
        self.description = str(description)
        self.product = product

    def __str__(self) -> str:
        if not self.product:
            return "%s (%s) %s" % (self.function, self.tool, self.description)
        return "%s (%s) %s: %s" % (self.function, self.tool, self.product, self.description)

    def __repr__(self) -> str:
        fmt = "GeneFunctionAnnotation(function=%r, tool='%s', product='%s', '%s')"
        return fmt % (self.function, self.tool, self.product, self.description)

    def __eq__(self, other: Any) -> bool:
        return (isinstance(other, _GeneFunctionAnnotation)
                and self.function == other.function
                and self.tool == other.tool
                and self.description == other.description
                and self.product == other.product)

    @classmethod
    def from_string(cls, text: str) -> "_GeneFunctionAnnotation":
        """ Regenerates a GeneFunctionAnnotation from a string version """
        try:
            parts = _parse_format("{} ({}) {}: {}", text)
        except ValueError:
            try:
                parts = _parse_format("{} ({}) {}", text)
            except ValueError:
                raise ValueError("cannot parse GeneFunctionAnnotation from %s" % text)

        if len(parts) not in [3, 4]:
            raise ValueError("cannot parse GeneFunctionAnnotation from %s" % text)

        if len(parts) == 4:
            function, tool, product, description = parts
            return cls(GeneFunction.from_string(function), tool, description, product)

        function, tool, description = parts
        return cls(GeneFunction.from_string(function), tool, description, product=None)


class GeneFunctionAnnotations:
    """ Tracks multiple annotations of gene functions. Each annotation contains
        the declared function of the gene, which tool determined the function,
        and a comment on how the determination was made.
    """
    slots = ["_annotations", "_by_tool", "_by_function"]

    def __init__(self) -> None:
        self._annotations = []  # type: List["_GeneFunctionAnnotation"]
        self._by_tool = defaultdict(list)  # type: Dict[str, List["_GeneFunctionAnnotation"]]
        self._by_function = defaultdict(list)  # type: Dict[GeneFunction, List["_GeneFunctionAnnotation"]]

    def __iter__(self) -> Iterator["_GeneFunctionAnnotation"]:
        for annotation in self._annotations:
            yield annotation

    def __len__(self) -> int:
        return len(self._annotations)

    def add(self, function: GeneFunction, tool: str, description: str,
            product: Optional[str] = None) -> "_GeneFunctionAnnotation":
        """ Adds a gene function annotation. If an existing annotation has all
            the same values as those provided, a duplicate will not be added.

            Arguments:
                function: a GeneFunction value
                tool: a string naming the tool that determined the gene function
                description: description text for storing extra information
                product: the name of the product the annotation is related to,
                         required only if the function is CORE

            Returns:
                the GeneFunction added (or the existing one if duplicated)
        """
        tool = str(tool)
        description = str(description)
        # if there's already an exactly similar function annotation, skip adding
        new = _GeneFunctionAnnotation(function, tool, description, product)
        existing_functions = self._by_function.get(function, [])
        for existing in existing_functions:
            if existing == new:
                return existing
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
            annotation = _GeneFunctionAnnotation.from_string(section)
            self.add(annotation.function, annotation.tool, annotation.description, annotation.product)

    def get_by_tool(self, tool: str) -> List:
        """ Returns a list of all GeneFunctionAnnotations which were added with
            the tool name provided.
        """
        return self._by_tool.get(tool, [])

    def get_by_function(self, function: GeneFunction) -> List[_GeneFunctionAnnotation]:
        """ Returns a list of GeneFunctionAnnotations which have the same
            gene function as that provided.
        """
        return self._by_function.get(function, [])

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
        # then priority for resfam, then smcogs
        annotations = self._by_tool.get("resfam", self._by_tool.get("smcogs"))
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
