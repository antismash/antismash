# License: GNU Affero General Public License v3 or later
# A copy of GNU AGPL v3 should have been included in this software package in LICENSE.txt.

""" A collection of classes and functions for handling rule categories
"""

from dataclasses import dataclass
from typing import Any, IO, Optional, Type, TypeVar, Union

from antismash.common import json


TRuleCategory = TypeVar("TRuleCategory", bound="RuleCategory")  # pylint: disable=invalid-name


@dataclass
class RuleCategory:
    """A collection of information about a particular rule category"""
    name: str
    description: str
    version: int = 1
    parent: Optional["RuleCategory"] = None

    @classmethod
    def from_json(cls: Type[TRuleCategory], name: str,
                  entry: dict[str, Union[str, int]],
                  existing_categories: dict[str, TRuleCategory]) -> TRuleCategory:
        """Regenerates a RuleCategory instance from a JSON representation"""

        try:
            description: str = str(entry["description"])  # cast to make mypy happy
            version: int = int(entry["version"])  # cast to make mypy happy
        except KeyError as err:
            raise ValueError(f"missing required rule metadata key '{err}'")

        parent: Optional[RuleCategory] = None
        if "parent" in entry:
            parent = existing_categories[str(entry["parent"])]

        if version != 1:
            raise ValueError(f"Unknown rule category version {version}")

        return cls(name, description, version, parent=parent)


def _parse_categories(data: dict[str, Any]) -> list[RuleCategory]:
    """ Parses categories out from a JSON-like dictionary structure

        Arguments:
            data: the raw JSON data

        Returns:
            a dictionary mapping category name to RuleCategory instance
    """
    # split entries into two sets, those without parents, then those with parents
    top_level = {}
    children = {}
    for name, metadata in data.items():
        if "parent" in metadata:
            if metadata["parent"] not in data:
                raise ValueError(f"category {name!r} refers to unknown category {metadata['parent']!r}")
            children[name] = metadata
        else:
            top_level[name] = metadata

    categories: dict[str, RuleCategory] = {}

    # run on the top level first, to ensure that they exist to be referenced
    for group in [top_level, children]:
        for name, metadata in group.items():
            categories[name] = RuleCategory.from_json(name, metadata, existing_categories=categories)

    return list(categories.values())


_PARSED_FILES: dict[str, list[RuleCategory]] = {}


def parse_category_file(handle: Union[str, IO]) -> list[RuleCategory]:
    """ Parses categories out from a JSON-like dictionary structure

        Arguments:
            data: the raw JSON data

        Returns:
            a dictionary mapping category name to RuleCategory instance
    """
    # if called before, just return the cached data
    if isinstance(handle, str):
        filename = handle
    else:
        filename = handle.name

    if filename in _PARSED_FILES:
        return _PARSED_FILES[filename]

    # not cached, generate
    if isinstance(handle, str):
        with open(filename, "r", encoding="utf-8") as real_handle:
            data = json.load(real_handle)
    else:
        data = json.load(handle)

    categories = _parse_categories(data)
    # and cache for future calls
    _PARSED_FILES[filename] = categories

    return categories
