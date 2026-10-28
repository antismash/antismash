# License: GNU Affero General Public License v3 or later
# A copy of GNU AGPL v3 should have been included in this software package in LICENSE.txt.

""" Ruleset construction for the subcluster detection module """

from typing import cast

from antismash.common import path
from antismash.common.hmm_rule_parser.cluster_prediction import Ruleset, create_rules

from .signatures import AGGREGATE_HMM_FILE, HmmSignature, get_signatures

# Strictness levels are cumulative: each level includes all levels above it,
# ordered from most to least strict.
STRICTNESS_LEVELS = ("strict", "relaxed", "loose")

# Categories referenced across all rule files, should always be "Subcluster"
_CATEGORIES = {"Subcluster"}

_CACHE: dict[str, Ruleset] = {}


def get_ruleset(strictness: str) -> Ruleset:
    """ Returns the parsed rules and signatures for the given strictness level

        Results are cached, with a single ruleset built per strictness level.

        Arguments:
            strictness: one of the supported strictness levels, from most to
                        least strict

        Returns:
            a ruleset ready for use by the rule-based detection pipeline
    """
    if strictness not in STRICTNESS_LEVELS:
        raise ValueError(f"Unknown strictness level {strictness!r}. "
                         f"Valid levels: {', '.join(STRICTNESS_LEVELS)}")

    if strictness in _CACHE:
        return _CACHE[strictness]

    ruleset = _build_ruleset(strictness)
    _CACHE[strictness] = ruleset
    return ruleset


def _rule_files_for_strictness(strictness: str) -> list[str]:
    """ Returns the ordered list of rule files to load for a given strictness level

        Files are ordered from the most strict level to the requested level, so
        that rules from stricter levels are always included.

        Arguments:
            strictness: the strictness level to gather rule files for

        Returns:
            a list of paths to rule files
    """
    levels = STRICTNESS_LEVELS[: STRICTNESS_LEVELS.index(strictness) + 1]
    return [path.get_full_path(__file__, "subcluster_rules", f"{level}.txt") for level in levels]


def _build_ruleset(strictness: str) -> Ruleset:
    """ Builds, without caching, the ruleset for the given strictness level

        Arguments:
            strictness: the strictness level to build the ruleset for

        Returns:
            a ruleset ready for use by the rule-based detection pipeline
    """
    signatures = get_signatures()
    rules = create_rules(_rule_files_for_strictness(strictness), set(signatures), _CATEGORIES)

    return Ruleset(
        tuple(rules),
        cast(dict[str, HmmSignature], signatures),
        AGGREGATE_HMM_FILE,
        _CATEGORIES,
        "subclusters",
        equivalence_groups=[],
    )
