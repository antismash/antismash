# License: GNU Affero General Public License v3 or later
# A copy of GNU AGPL v3 should have been included in this software package in LICENSE.txt.

""" JSON-friendly classes explicitly for use by the javascript drawing libraries
"""

from dataclasses import dataclass, field
from json import JSONDecodeError  # pylint: disable=unused-import  # for compatibility
from typing import Any, Callable, IO, Optional, Self, TypeAlias, Union

# pylint doesn't recognise any members of orjson, where mypy and pyflakes do
# pylint: disable=no-name-in-module
from orjson import (
    loads,  # pylint: disable=unused-import  # used by others
    OPT_NON_STR_KEYS as _OPT_CONVERT_NON_STR_KEYS,
    OPT_SORT_KEYS as _OPT_SORT_KEYS,
    OPT_INDENT_2 as _OPT_INDENT_2,
    dumps as _dumps,
)
# pylint: enable=no-name-in-module


from antismash.common.secmet.record import Seq
from antismash.common.secmet.features import CDSFeature
from antismash.common.secmet.qualifiers import NRPSPKSQualifier


@dataclass(slots=True)
class JSONBase:
    """ A base class for JSON-serialisable objects """


JSONCompatible: TypeAlias = Union[
    dict[str, Union["JSONCompatible"]],
    list["JSONCompatible"],
    tuple["JSONCompatible"],
    JSONBase,
    int,
    float,
    str,
    bool,
    None
]


def _base_convertor(obj: Any) -> Any:
    # handles any conversion methods for classes that aren't default types or dataclasses
    if isinstance(obj, Seq):
        return str(obj)
    if hasattr(obj, "to_json"):
        return obj.to_json()
    if hasattr(obj, "__json__"):
        return obj.__json__()
    # but if no conversion method is found, then an error must be raised for orjson (and stdlib json, for that matter)
    raise TypeError


def _convert_std_to_orson(*, sort_keys: bool = False, option: int = 0, indent: bool = True) -> int:
    # always match stdlib JSON's default behaviour, where non-string keys are converted to string
    option |= _OPT_CONVERT_NON_STR_KEYS
    if sort_keys:
        option |= _OPT_SORT_KEYS
    if indent:
        option |= _OPT_INDENT_2
    return option


def dump(obj: Any, handle: IO, *, default: Callable[[Any], Any] = _base_convertor, indent: bool = False,
         sort_keys: bool = False, option: int = 0,
         ) -> None:
    """ Converts the given object to JSON and writes the resulting string to the given file handle

        Arguments:
            obj: the object to convert
            handle: the file object to write to
            default: an optional override of the usual class convertor handler for non-standard types
            indent: a boolean indicating whether to use indents in the string conversion (always 2 spaces if used)
            sort_keys: whether the child attributes should be sorted by key
            option: an orjson option value (see orjson documentation for possible values)

        Returns:
            None
    """
    option = _convert_std_to_orson(indent=indent, sort_keys=sort_keys, option=option)
    handle.write(dumps(obj, default=default, option=option))


def dumps(obj: Any, *, default: Callable[[Any], Any] = _base_convertor, indent: bool = False,
          sort_keys: bool = False, option: int = 0,
          ) -> str:
    """ Converts the given object to a JSON string

        Arguments:
            obj: the object to convert
            default: an optional override of the usual class convertor handler for non-standard types
            indent: a boolean indicating whether to use indents in the string conversion (always 2 spaces if used)
            sort_keys: whether the child attributes should be sorted by key
            option: an orjson option value (see orjson documentation for possible values)

        Returns:
            the string generated
    """
    option = _convert_std_to_orson(indent=indent, sort_keys=sort_keys, option=option)
    return _dumps(obj, default=default, option=option).decode()


def load(handle: IO) -> dict[str, Any]:
    """ Reads in JSON text from the given file handle and returns the information using
        standard types.

        Arguments:
            handle: the file handle to read from

        Returns:
            a dictionary mapping loaded key-value pairs
    """
    return loads(handle.read())


@dataclass
class JSONDomain(JSONBase):
    """ A JSON-serialisable object for simplifying domain datatypes throughout this file """
    type: str
    start: int
    end: int
    predictions: list[tuple[str, str]]
    napdoslink: str
    blastlink: str
    sequence: str
    dna: str
    abbreviation: str
    html_class: str

    @classmethod
    def from_domain(cls, domain: NRPSPKSQualifier.Domain, *args: Any) -> Self:
        """ A helper for constructing an instance from an existing domain """
        return cls(domain.full_type, domain.start, domain.end, *args)


@dataclass
class JSONModule(JSONBase):
    """ A JSON-serialisable object for simplifying NRPS/PKS module datatypes """
    start: int
    end: int
    complete: bool
    iterative: bool
    monomer: str
    multi_cds: Optional[str] = None
    match_id: Optional[str] = None

    def __post_init__(self) -> None:
        if bool(self.multi_cds) != bool(self.match_id):
            raise ValueError("multi_cds and match_id must both have values or both be None")


@dataclass
class JSONOrf(JSONBase):
    """ A JSON-serialisable object for simplifying ORF datatypes throughout this file """
    sequence: str
    id: str
    domains: list[JSONDomain] = field(default_factory=list)
    modules: list[JSONModule] = field(default_factory=list)

    @classmethod
    def from_cds(cls, feature: CDSFeature) -> Self:
        """ A helper for constructing an instance from an existing CDS feature """
        return cls(feature.translation, feature.get_name())

    def add_domain(self, domain: JSONDomain) -> None:
        """ Add a JSONDomain to the list of domains in this ORF """
        assert isinstance(domain, JSONDomain)
        self.domains.append(domain)

    def add_module(self, module: JSONModule) -> None:
        """ Add a JSONModule to the list of modules contained by this ORF """
        assert isinstance(module, JSONModule)
        self.modules.append(module)
