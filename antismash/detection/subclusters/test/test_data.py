# License: GNU Affero General Public License v3 or later
# A copy of GNU AGPL v3 should have been included in this software package in LICENSE.txt.

# for test files, silence irrelevant and noisy pylint warnings
# pylint: disable=use-implicit-booleaness-not-comparison,protected-access,missing-docstring

import os
import unittest

from antismash.common import path
from antismash.common.hmm_rule_parser.rule_parser import Tokeniser, TokenTypes
from antismash.detection.subclusters import compounds
from antismash.detection.subclusters.ruleset import STRICTNESS_LEVELS, get_ruleset
from antismash.detection.subclusters.signatures import get_signatures

COMPOUND_FILE = path.get_full_path(compounds.__file__, "data", "compound_details.txt")


class TestDataConsistency(unittest.TestCase):
    def setUp(self):
        self.compounds = compounds._read_compounds(COMPOUND_FILE)

    def test_every_rule_has_compound(self):
        for strictness in STRICTNESS_LEVELS:
            missing = get_ruleset(strictness).get_rule_names() - set(self.compounds)
            assert not missing, f"rules without compound in {strictness}: {sorted(missing)}"

    def test_every_compound_has_rule(self):
        # the least strict level contains the rules of all levels
        rule_names = get_ruleset(STRICTNESS_LEVELS[-1]).get_rule_names()
        unused = set(self.compounds) - rule_names
        assert not unused, f"compounds without rule: {sorted(unused)}"

    def test_signatures_match_hmm_files(self):
        for name, signature in get_signatures().items():
            assert os.path.exists(signature.hmm_file), signature.hmm_file
            with open(signature.hmm_file, encoding="utf-8") as handle:
                hmm_names = [line.split()[-1] for line in handle if line.startswith("NAME")]
            assert hmm_names == [name], f"{signature.hmm_file}: {hmm_names} != [{name!r}]"

    def test_extenders_have_signatures(self):
        # the rule parser only checks that the identifiers in a rule's conditions
        # have signatures, so an extender without one would silently never match
        signatures = set(get_signatures())
        missing = {}
        for rule in get_ruleset(STRICTNESS_LEVELS[-1]).rules:
            if not rule.extenders:
                continue
            tokens = Tokeniser(str(rule.extenders)).tokens
            names = {token.identifier for token in tokens if token.type == TokenTypes.IDENTIFIER}
            if names - signatures:
                missing[rule.name] = sorted(names - signatures)
        assert not missing, f"extenders without signatures: {missing}"
