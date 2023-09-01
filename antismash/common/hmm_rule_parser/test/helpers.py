# License: GNU Affero General Public License v3 or later
# A copy of GNU AGPL v3 should have been included in this software package in LICENSE.txt.

# for test files, silence irrelevant and noisy pylint warnings
# pylint: disable=protected-access,missing-docstring,too-many-public-methods

""" helper functions for testing modules using hmm_rule_parser """

import os

from antismash.common.hmm_rule_parser.cluster_prediction import Ruleset
from antismash.common.signature import get_signature_profiles


def check_hmm_signatures(signature_file, hmm_dir):
    sigs = get_signature_profiles(signature_file)
    for sig in sigs:
        hmm_file = os.path.abspath(os.path.join(hmm_dir, sig.path))
        assert os.path.exists(hmm_file)
        name = None
        with open(hmm_file, encoding="utf-8") as handle:
            for line in handle:
                if line.startswith("NAME"):
                    name = line.split()[-1]
        assert name
        assert name == sig.name, f"{name} != {sig.name}"


def create_ruleset(rules, *, categories=None, hmm_profiles=None, seeds="dummy_seeds", tool="test_tool",
                   dynamic_profiles=None, equivalence_groups=None):
    if categories is None:
        categories = {rule.category for rule in rules}
    return Ruleset(tuple(rules), hmm_profiles or {}, seeds, categories,
                   tool=tool, dynamic_profiles=dynamic_profiles or {},
                   equivalence_groups=equivalence_groups or set(),
                   )
