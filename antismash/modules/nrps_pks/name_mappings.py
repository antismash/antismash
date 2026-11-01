# License: GNU Affero General Public License v3 or later
# A copy of GNU AGPL v3 should have been included in this software package in LICENSE.txt.
# pylint: disable=line-too-long
"""A-domain substrate names, structures, and taxonomy parent relationships.
"""

from __future__ import annotations

from dataclasses import dataclass, field, fields
from typing import Any

from antismash.common import json, path

CATEGORY_FILE = path.get_full_path(__file__, "data", "substrate_categories.json")
SUBSTRATE_FILE = path.get_full_path(__file__, "data", "substrates.json")


def _ancestors(parent: Category | None) -> tuple[Category, ...]:
    result = []
    while parent is not None:
        result.append(parent)
        parent = parent.parent
    return tuple(result)


@dataclass(frozen=True)
class Category:
    """A taxonomy category linked to its parent, ending at a gross-class root.

    smiles is an illustrative common scaffold with wildcard attachment sites;
    these sites can be hydrogen in some members. Use consensus_smarts to query
    the conserved core, and parent links to determine taxonomy membership.
    """

    id: str
    name: str
    level: str
    parent: Category | None = field(default=None, repr=False)
    smiles: str | None = None
    consensus_smarts: str | None = None
    consensus_status: str = "unavailable"

    @property
    def ancestors(self) -> tuple[Category, ...]:
        """Ancestors from the immediate parent to the gross-class root."""
        return _ancestors(self.parent)

    @property
    def children(self) -> tuple[Category | Substrate, ...]:
        """Immediate child categories or substrates, in data-file order."""
        return _CHILDREN_BY_PARENT.get(self.id, ())

    @property
    def substrates(self) -> tuple[Substrate, ...]:
        """All descendant substrates, traversing children in data-file order."""
        return tuple(
            leaf for child in self.children
            for leaf in (child.substrates if isinstance(child, Category) else (child,))
        )

    @classmethod
    def from_json(cls, existing: dict[str, Category], parent_id: str | None = None, **kwargs: Any) -> Category:
        """ Load from JSON data, linking to an existing parent if provided. """
        parent = None
        if parent_id is not None:
            parent = existing[parent_id]
        return cls(parent=parent, **kwargs)

    def to_json(self) -> dict[str, Any]:
        """Serialize parent links as parent_id, without recursively duplicating parents."""
        return {
            ("parent_id" if item.name == "parent" else item.name): (self.parent.id if self.parent else None)
            if item.name == "parent" else getattr(self, item.name)
            for item in fields(self)
        }

    def __str__(self) -> str:
        return self.name


def load_categories() -> dict[str, Category]:
    """ Load substrate categories from the JSON file and build the parent-child relationships. """
    categories_by_id: dict[str, Category] = {}

    with open(CATEGORY_FILE, "r", encoding="utf-8") as handle:
        category_data: list[dict[str, str | None]] = json.load(handle)["data"]

    for cat in category_data:
        category = Category.from_json(categories_by_id, **cat)
        categories_by_id[category.id] = category

    return categories_by_id


CATEGORIES_BY_ID: dict[str, Category] = load_categories()


@dataclass(frozen=True, order=True)
class Substrate:  # pylint: disable=too-many-instance-attributes
    """Names and taxonomy-sourced structures for a terminal substrate.

    Structurally unresolved entries can contain wildcard atoms.
    A None stereochemical_smiles means no resolved stereochemical structure was
    recorded; it is never replaced by the connectivity-only SMILES.
    """

    long: str
    short: str
    norine: str
    substrate_type: str
    smiles: str | None = field(default=None, compare=False)
    stereochemical_smiles: str | None = field(default=None, compare=False)
    id: str | None = field(default=None, compare=False)
    parent: Category | None = field(default=None, repr=False, compare=False)

    @property
    def ancestors(self) -> tuple[Category, ...]:
        """Ancestors from the small cluster through to the gross class."""
        return _ancestors(self.parent)

    def to_json(self) -> dict[str, Any]:
        """Serialize all fields, representing the parent by category ID."""
        return {
            item.name: (self.parent.id if self.parent else None)
            if item.name == "parent" else getattr(self, item.name)
            for item in fields(self)
        }

    @classmethod
    def from_json(cls, data: dict[str, Any]) -> Substrate:
        """Load names and metadata, accepting file parent_id and serialized parent links."""
        values = {item.name: data[item.name] for item in fields(cls) if item.name in data}
        parent_id = values.pop("parent", data.get("parent_id"))
        if parent_id is not None:
            values["parent"] = CATEGORIES_BY_ID[parent_id]
        return cls(**values)

    def __str__(self) -> str:
        """Use the antiSMASH short name as the canonical string representation."""
        return self.short


KNOWN_CATEGORIES: tuple[Category, ...] = tuple(CATEGORIES_BY_ID.values())


def load_substrates() -> list[Substrate]:
    """Load substrate names, structures, and category links from the JSON file."""
    with open(SUBSTRATE_FILE, "r", encoding="utf-8") as handle:
        substrate_data = json.load(handle)["data"]

    return [Substrate.from_json(record) for record in substrate_data]


KNOWN_SUBSTRATES: list[Substrate] = load_substrates()

SUBSTRATES_BY_ID: dict[str, Substrate] = {
    substrate.id: substrate for substrate in KNOWN_SUBSTRATES if substrate.id is not None
}
_children: dict[str, list[Category | Substrate]] = {}
for _category in KNOWN_CATEGORIES:
    if _category.parent is not None:
        _children.setdefault(_category.parent.id, []).append(_category)
for _substrate in KNOWN_SUBSTRATES:
    if _substrate.parent is not None:
        _children.setdefault(_substrate.parent.id, []).append(_substrate)
_CHILDREN_BY_PARENT = {key: tuple(value) for key, value in _children.items()}

_LONG_TO_SUBSTRATE: dict[str, Substrate] = {
    sub.long.lower(): sub for sub in KNOWN_SUBSTRATES
}
_SHORT_TO_SUBSTRATE: dict[str, Substrate] = {
    sub.short.lower(): sub for sub in KNOWN_SUBSTRATES
}
_NORINE_TO_SUBSTRATE: dict[str, Substrate] = {
    sub.norine.lower(): sub for sub in KNOWN_SUBSTRATES
}


def get_substrate_by_name(name: str) -> Substrate:
    """Look up a substrate by taxonomy ID, long, short, or Norine name.
    """
    if name in SUBSTRATES_BY_ID:
        return SUBSTRATES_BY_ID[name]
    lc_name = name.lower()
    for index in (_LONG_TO_SUBSTRATE, _SHORT_TO_SUBSTRATE, _NORINE_TO_SUBSTRATE):
        if lc_name in index:
            return index[lc_name]
    raise ValueError(f"Substrate {name!r} not found")


def get_category_by_name(name: str) -> Category:
    """Look up an exact model/taxonomy ID or a case-insensitive display name.
    """
    if name in CATEGORIES_BY_ID:
        return CATEGORIES_BY_ID[name]
    matches = [category for category in KNOWN_CATEGORIES if category.name.lower() == name.lower()]
    if len(matches) == 1:
        return matches[0]
    if matches:
        raise ValueError(f"Ambiguous category {name!r}; use a taxonomy ID")
    raise ValueError(f"Category {name!r} not found")


def lookup(name: str) -> Substrate | Category:
    """Look up a substrate or category by name or ID."""
    try:
        return get_substrate_by_name(name)
    except ValueError:
        return get_category_by_name(name)


def is_valid_norine_name(name: str) -> bool:
    """ Check if a name is valid for Norine """
    return name.lower() in _NORINE_TO_SUBSTRATE
