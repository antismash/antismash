# License: GNU Affero General Public License v3 or later
# A copy of GNU AGPL v3 should have been included in this software package in LICENSE.txt.

""" A thin wrapper around antismash.common.hmm_rule_parser.categories to avoid
    needing path information
"""

from antismash.common.path import get_full_path
from antismash.common.hmm_rule_parser.categories import parse_category_file, RuleCategory

_CATEGORIES: list[RuleCategory] = []


def get_rule_categories() -> list[RuleCategory]:
    """Generate a list of rule categories from categories.json
    Only processes the file once per Python invocation, future calls access the cached data.
    """
    if not _CATEGORIES:
        _CATEGORIES.extend(parse_category_file(get_full_path(__file__, "data", "categories.json")))
    return _CATEGORIES
