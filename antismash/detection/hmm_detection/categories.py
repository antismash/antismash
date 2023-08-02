# License: GNU Affero General Public License v3 or later
# A copy of GNU AGPL v3 should have been included in this software package in LICENSE.txt.

""" Handle valid rule categories plus associated metadata
"""

from dataclasses import dataclass
import json
from typing import Any, Dict, List, Type, Optional, TypeVar, Union

from antismash.common import path

TRuleCategory = TypeVar("TRuleCategory", bound="RuleCategory")


@dataclass
class RuleCategory:
    """A collection of information about a particular rule category"""
    name: str
    description: str
    version: int = 1
    parent: Optional["RuleCategory"] = None

    @classmethod
    def from_json(cls: Type[TRuleCategory], name: str,
                  entry: Dict[str, Union[str, int]],
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


def _parse_categories(data: dict[str, Any]) -> dict[str, RuleCategory]:
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

    return categories


def get_rule_categories() -> List[RuleCategory]:
    """Generate a list of rule categories from categories.json
    Only processes the file once per Python invocation, future calls access the cached data.
    """
    # if called before, just return the cached data
    existing = getattr(get_rule_categories, 'existing', None)
    if existing is not None:
        assert isinstance(existing, list)
        return existing

    # not cached, generate
    category_mapping: dict[str, RuleCategory] = {}

    with open(path.get_full_path(__file__, "data", "categories.json"), "r", encoding="utf-8") as handle:
        category_mapping = _parse_categories(json.load(handle))

    categories = list(category_mapping.values())
    # and cache for future calls, and silence mypy warning as mypy can't handle this
    get_rule_categories.existing = categories  # type: ignore

    return categories
