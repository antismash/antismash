# License: GNU Affero General Public License v3 or later
# A copy of GNU AGPL v3 should have been included in this software package in LICENSE.txt.

""" Data handling functions and classes for the terpene module
"""
from dataclasses import asdict, dataclass
from typing import Any, Optional, Self

from antismash.common import json, path


class MissingCompoundError(ValueError):
    """ A specific error for when expected compounds are missing """


class MissingHmmError(ValueError):
    """ A specific error for when expected profiles are missing """


@dataclass(frozen=True, slots=True)
class CompoundGroup:
    """ Biosynthetic and chemical properties for a group of compounds.
    """
    name: str
    extended_name: Optional[str]
    single_compound: bool
    biosynthetic_class: str
    chain_length: int
    initial_cyclisations: tuple[str]
    functional_groups: tuple[str]
    biosynthetic_subclass: str = ""

    def get_compound_name(self) -> str:
        """ Returns the name of the compound, if the compound group
            represents a single compound, otherwise returns an empty string
        """
        if not self.single_compound:
            return ""
        return self.name

    def get_cyclisations_description(self) -> str:
        """ Returns a description of the cyclisations for this
            compound group
        """
        if not self.initial_cyclisations:
            return "acyclic"
        if "unknown" in self.initial_cyclisations:
            return "unknown"
        return " + ".join(self.initial_cyclisations)

    def get_functional_groups_description(self) -> str:
        """ Returns the functional groups, if it exists,
            otherwise returns "none"
        """
        if not self.functional_groups:
            return "none"
        if "unknown" in self.functional_groups:
            return "unknown"
        return ", ".join(self.functional_groups)

    def to_json(self) -> dict[str, Any]:
        """ Returns a JSON-friendly representation """
        return asdict(self)

    def __hash__(self) -> int:
        return hash(self.name)

    @classmethod
    def from_json(cls, data: dict[str, Any]) -> Self:
        """ Reconstructs an instance from a JSON representation """
        data["initial_cyclisations"] = tuple(data["initial_cyclisations"])
        data["functional_groups"] = tuple(data["functional_groups"])
        data["biosynthetic_subclass"] = data["biosynthetic_subclass"] or ""
        return cls(**data)


class Reaction:
    """ Contains the substrates and products of a chemical reaction.
    """
    def __init__(self, substrates: tuple[CompoundGroup, ...], products: tuple[CompoundGroup, ...]):
        self.substrates = substrates
        self.products = products

    def __eq__(self, other: Any) -> bool:
        if not isinstance(other, Reaction):
            return False
        return self.substrates == other.substrates and self.products == other.products

    def build_intersection(self, other: "Reaction") -> "Reaction":
        """ Creates a new Reaction instance that contains the
            intersections of the substrates and products of this instance
            and the provided instance.
        """
        return Reaction(
            tuple(set(self.substrates) & set(other.substrates)),
            tuple(set(self.products) & set(other.products)),
        )

    def __str__(self) -> str:
        return (f"substrates={[compound.name for compound in self.substrates]}, "
                f"products={[compound.name for compound in self.products]}")

    def to_json(self) -> dict[str, Any]:
        """ Returns a JSON-friendly representation  """
        return {"substrates": [compound.name for compound in self.substrates],
                "products": [compound.name for compound in self.products]}

    @classmethod
    def from_json(cls, data: dict[str, list[str]], compound_groups: dict[str, CompoundGroup],
                  ) -> Self:
        """ Reconstructs an instance from a JSON representation """
        missing: list[str] = []
        for key in ["substrates", "products"]:
            missing.extend(item for item in data[key] if item not in compound_groups)
        if missing:
            raise MissingCompoundError(f"Compound groups not defined: {', '.join(missing)}")

        return cls(tuple(compound_groups[name] for name in data["substrates"]),
                   tuple(compound_groups[name] for name in data["products"]))


@dataclass
class TerpeneHMM:
    """ Properties associated with a terpene hmm profile
    """
    name: str
    description: str
    domain_type: str
    length: int
    cutoff: int
    subtypes: tuple["TerpeneHMM", ...]
    reactions: tuple[Reaction, ...]
    __is_subtype: bool = False

    def __post_init__(self) -> None:
        for subtype in self.subtypes:
            subtype.mark_as_subtype()
        if self.reactions:
            first = self.reactions[0]
            for reaction in self.reactions[1:]:
                if any(substrate in first.substrates for substrate in reaction.substrates):
                    raise ValueError("All profile reactions must have mutually exclusive substrates")

    def is_subtype(self) -> bool:
        """ Returns whether the instance is a subtype """
        return self.__is_subtype

    def mark_as_subtype(self) -> None:
        """ Marks the instance as a subtype """
        self.__is_subtype = True

    @classmethod
    def from_json(cls, hmm_json: dict[str, Any], terpene_hmms: dict[str, "TerpeneHMM"],
                  compound_groups: dict[str, CompoundGroup]) -> Self:
        """ Reconstructs an instance from a JSON representation """
        try:
            subtypes = tuple(terpene_hmms[name] for name in hmm_json["subtypes"])
        except KeyError as key:
            raise MissingHmmError(f"'{hmm_json['name']}': Subtype {key} not defined yet")
        return cls(str(hmm_json["name"]), str(hmm_json["description"]), str(hmm_json["type"]),
                   int(hmm_json["length"]), int(hmm_json["cutoff"]), subtypes,
                   tuple(Reaction.from_json(reaction, compound_groups)
                         for reaction in hmm_json["reactions"]))


_COMPOUND_CACHE: dict[str, CompoundGroup] = {}
_HMM_PROPERTIES_CACHE: dict[str, TerpeneHMM] = {}
_HMM_LENGTHS: dict[str, int] = {}
HMM_METADATA_FILE = path.get_full_path(__file__, "data", "hmm_properties.json")
COMPOUND_FILE = path.get_full_path(__file__, "data", "compound_groups.json")


def load_compounds() -> dict[str, CompoundGroup]:
    """ Loads compound groups from compound_groups.json
        Only does the processing once per python invocation, future runs access
        existing compounds

        Arguments:
            None

        Returns:
            a dictionary of compound group names to compound groups
    """
    if _COMPOUND_CACHE:
        return _COMPOUND_CACHE
    compound_groups: dict[str, CompoundGroup] = {}
    with open(COMPOUND_FILE, encoding="utf-8") as handle:
        compounds_json = json.load(handle)
    for group_data in compounds_json["groups"]:
        compound_group = CompoundGroup.from_json(group_data)
        compound_groups[compound_group.name] = compound_group
    _COMPOUND_CACHE.update(compound_groups)
    return _COMPOUND_CACHE


def load_hmm_properties() -> dict[str, TerpeneHMM]:
    """ Loads the HMM properties from hmm_properties.json
        Only does the processing once per python invocation, future runs access
        existing properties

        Arguments:
            None

        Returns:
            a dictionary of hmm names to TerpeneHMM objects
    """
    # if already called once, then just reuse the cached results
    if not _HMM_PROPERTIES_CACHE:
        with open(HMM_METADATA_FILE, encoding="utf-8") as handle:
            properties_json = json.load(handle)

        terpene_hmms: dict[str, TerpeneHMM] = {}
        for hmm_data in properties_json["profiles"]:
            terpene_hmm = TerpeneHMM.from_json(hmm_data, terpene_hmms, load_compounds())
            terpene_hmms[terpene_hmm.name] = terpene_hmm

        _HMM_PROPERTIES_CACHE.update(terpene_hmms)

    return _HMM_PROPERTIES_CACHE


def load_hmm_lengths(hmm_properties: dict[str, TerpeneHMM]) -> dict[str, int]:
    """ Loads hmm lengths from hmm properties.
        Only does the processing once per python invocation, future runs access
        existing properties

        Arguments:
            hmm_properties: a dictionary of hmm names to TerpeneHMM objects

        Returns:
            a dictionary of hmm names to hmm lengths
    """
    if not _HMM_LENGTHS:
        hmm_lengths = {hmm_name: hmm_obj.length for hmm_name, hmm_obj in hmm_properties.items()}
        _HMM_LENGTHS.update(hmm_lengths)
    return _HMM_LENGTHS
