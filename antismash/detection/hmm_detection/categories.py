# License: GNU Affero General Public License v3 or later
# A copy of GNU AGPL v3 should have been included in this software package in LICENSE.txt.

""" Handle valid rule categories plus associated metadata
"""

from dataclasses import dataclass
import json
from typing import Dict, List, Type, TypeVar, Union

from antismash.common import path

TRuleCategory = TypeVar("TRuleCategory", bound="RuleCategory")


@dataclass
class RuleCategory:
    """A collection of information about a particular rule category"""
    name: str
    description: str
    version: int = 1

    @classmethod
    def from_json(cls: Type[TRuleCategory], name: str,
                  entry: Dict[str, Union[str, int]]) -> TRuleCategory:
        """Regenerates a RuleCategory instance from a JSON representation"""

        try:
            description: str = str(entry["description"])  # cast to make mypy happy
            version: int = int(entry["version"])  # cast to make mypy happy
        except KeyError as err:
            raise ValueError(f"missing required rule metadata key '{err}'")

        if version != 1:
            raise ValueError(f"Unknown rule category version {version}")

        return cls(name, description, version)


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
    categories: List[RuleCategory] = []

    with open(path.get_full_path(__file__, "data", "categories.json"), 'r') as handle:
        category_json = json.load(handle)
        for name, metadata in category_json.items():
            categories.append(RuleCategory.from_json(name, metadata))

    # and cache for future calls, and silence mypy warning as mypy can't handle this
    get_rule_categories.existing = categories  # type: ignore

    return categories
