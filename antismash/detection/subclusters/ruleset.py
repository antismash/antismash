# License: GNU Affero General Public License v3 or later
# A copy of GNU AGPL v3 should have been included in this software package in LICENSE.txt.

""" Ruleset construction for the subcluster detection module """

from typing import cast

from antismash.common import path
from antismash.common.hmm_rule_parser.cluster_prediction import Ruleset, create_rules

from .signatures import AGGREGATE_HMM_FILE, HmmSignature, get_signatures

RULES_FILE = path.get_full_path(__file__, "data", "subcluster_rules.txt")

# Categories referenced in the rule file, should always be "Subcluster"
_CATEGORIES = {"Subcluster"}

# keyed by rule file, so that a patched RULES_FILE gets its own ruleset
_CACHE: dict[str, Ruleset] = {}


def get_ruleset() -> Ruleset:
    """ Returns the parsed rules and signatures

        The ruleset is built once and cached.

        Returns:
            a ruleset ready for use by the rule-based detection pipeline
    """
    if RULES_FILE not in _CACHE:
        _CACHE[RULES_FILE] = _build_ruleset(RULES_FILE)
    return _CACHE[RULES_FILE]


def _build_ruleset(rule_file: str) -> Ruleset:
    """ Builds, without caching, the ruleset from the given rule file

        Arguments:
            rule_file: the path of the rule file to build the ruleset from

        Returns:
            a ruleset ready for use by the rule-based detection pipeline
    """
    signatures = get_signatures()
    rules = create_rules([rule_file], set(signatures), _CATEGORIES)

    return Ruleset(
        tuple(rules),
        cast(dict[str, HmmSignature], signatures),
        AGGREGATE_HMM_FILE,
        _CATEGORIES,
        "subclusters",
        equivalence_groups=[],
    )
