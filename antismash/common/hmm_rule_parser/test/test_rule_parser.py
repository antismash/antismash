# License: GNU Affero General Public License v3 or later
# A copy of GNU AGPL v3 should have been included in this software package in LICENSE.txt.

# for test files, silence irrelevant and noisy pylint warnings
# pylint: disable=no-self-use,protected-access,missing-docstring,too-many-public-methods

import unittest

from antismash.common.hmm_rule_parser import rule_parser
from antismash.common.test.helpers import DummyCDS, DummyRecord, FakeHSPHit


def format_as_rule(name, cutoff, extent, conditions):
    return "RULE {} CUTOFF {} EXTENT {} CONDITIONS {}".format(name, cutoff, extent, conditions)


class DetectionTest(unittest.TestCase):
    def setUp(self):
        self.feature_by_id = {
            "GENE_1": DummyCDS(0, 30000, locus_tag="GENE_1"),
            "GENE_2": DummyCDS(30000, 50000, locus_tag="GENE_2"),
            "GENE_3": DummyCDS(70000, 90000, locus_tag="GENE_3"),
            "GENE_X": DummyCDS(95000, 100000, locus_tag="GENE_X"),
            "GENE_4": DummyCDS(125000, 140000, locus_tag="GENE_4"),
            "GENE_5": DummyCDS(145000, 150000, locus_tag="GENE_5")
        }
        self.features = list(self.feature_by_id.values())
        self.features.sort(key=lambda x: x.location.start)  # vital for py3 < 3.5
        self.record = DummyRecord(self.features)

        self.results_by_id = {
            "GENE_1": [FakeHSPHit("a", "GENE_1", 0, 10, 50, 0),
                       FakeHSPHit("b", "GENE_1", 0, 10, 50, 0)],
            "GENE_2": [FakeHSPHit("a", "GENE_1", 0, 10, 50, 0),
                       FakeHSPHit("c", "GENE_1", 0, 10, 50, 0)],
            "GENE_3": [FakeHSPHit("b", "GENE_1", 0, 10, 50, 0),
                       FakeHSPHit("c", "GENE_1", 0, 10, 50, 0)],
            "GENE_4": [FakeHSPHit("e", "GENE_1", 0, 10, 50, 0),
                       FakeHSPHit("f", "GENE_1", 0, 10, 50, 0)],
            "GENE_5": [FakeHSPHit("f", "GENE_1", 0, 10, 50, 0),
                       FakeHSPHit("g", "GENE_1", 0, 10, 50, 0)]}
        self.signature_names = set(["a", "b", "c", "d", "e", "f", "g", "modelA", "modelB"])

    def run_test(self, name, cutoff, extent, conditions):
        rules = format_as_rule(name, cutoff, extent, conditions)
        rules = rule_parser.Parser(rules, self.signature_names).rules
        for rule in rules:
            assert rule.contains_positive_condition()

        detected_types = {}
        cds_with_hits = sorted(self.results_by_id,
                               key=lambda gene_id: self.feature_by_id[gene_id].location.start)
        for cds in cds_with_hits:
            assert cds not in detected_types
            rule_results = []
            for rule in rules:
                rule_results.append(rule.detect(cds, self.feature_by_id, self.results_by_id))
                if rule_results[-1]:
                    # check that we have something interesting to report
                    hit_string = rule.get_hit_string().replace("0*", "")
                    assert "*" in hit_string
            results = [rule.name for rule, res in zip(rules, rule_results) if res.met and res.matches]
            if results:
                detected_types[cds] = results

        for gid in detected_types:
            detected_types[gid] = set(detected_types[gid])
        return detected_types

    def expect(self, results, genes_to_hit):
        expected = {gene: set("A") for gene in genes_to_hit}
        assert results == expected

    def test_single(self):
        results = self.run_test("A", 10, 20, "a")
        self.expect(results, ["GENE_1", "GENE_2"])

    def test_simple_or(self):
        results = self.run_test("A", 10, 20, "a or c")
        self.expect(results, ["GENE_1", "GENE_2", "GENE_3"])

    def test_chained_or(self):
        results = self.run_test("A", 10, 20, "a or b or c or e or f")
        assert len(results) == 5  # all valid

    def test_simple_and(self):
        results = self.run_test("A", 10, 20, "a and c")
        self.expect(results, ["GENE_1", "GENE_2"])

        results = self.run_test("A", 10, 20, "e and c")
        self.expect(results, [])

    def test_cds_and(self):
        results = self.run_test("A", 10, 20, "cds(a and c)")
        # GENE_2 is an exact match
        # GENE_1 is in range, but doesn't satisfy any condition locally
        self.expect(results, ["GENE_2"])

    def test_cds_or(self):
        results = self.run_test("A", 10, 20, "cds(a or c)")
        self.expect(results, ["GENE_1", "GENE_2", "GENE_3"])

    def test_cds_internal_combo(self):
        results = self.run_test("A", 10, 20, "cds(a and (b or c))")
        self.expect(results, ["GENE_1", "GENE_2"])

    def test_chained_and_a(self):
        results = self.run_test("A", 10, 20, "a and b and not c")
        self.expect(results, ["GENE_1"])  # 2 reaches c

    def test_chained_and_b(self):
        results = self.run_test("A", 10, 20, "a and b and not cds(a and b)")
        # the range is too short to reach 2 for the a
        self.expect(results, [])

    def test_chained_and_c(self):
        results = self.run_test("A", 21, 20, "a and b and not cds(a and b)")
        self.expect(results, ["GENE_3"])  # 3 has b, 2 has a, 1 reaches 2 with a,b

    def test_simple_minimum(self):
        results = self.run_test("A", 10, 20, "minimum(2, [a, b])")
        self.expect(results, ["GENE_1", "GENE_2"])

        results = self.run_test("A", 100, 20, "minimum(2, [b, e])")
        self.expect(results, ["GENE_1", "GENE_3", "GENE_4"])

    def test_single_gene(self):
        self.results_by_id = {
            "GENE_1": [
                FakeHSPHit("modelA", "GENE_1", 0, 10, 50, 0),
                FakeHSPHit("modelB", "GENE_1", 0, 10, 50, 0)
            ]}
        self.feature_by_id = {
            "GENE_1": DummyCDS(0, 30000, locus_tag="GENE_1")
        }

        results = self.run_test("A", 10, 20, "minimum(2, [modelA,modelB])")
        self.expect(results, ["GENE_1"])

        results = self.run_test("A", 10, 20, "cds(modelA and modelB)")
        self.expect(results, ["GENE_1"])

        results = self.run_test("A", 10, 20, "cds(modelA or modelB)")
        self.expect(results, ["GENE_1"])

        results = self.run_test("A", 10, 20, "modelA and modelB")
        self.expect(results, ["GENE_1"])

        results = self.run_test("A", 10, 20, "modelA or modelB")
        self.expect(results, ["GENE_1"])

    def test_negated_minimum(self):
        # only GENE_2 and GENE_3 have a C
        results = self.run_test("A", 10, 20, "c and not minimum(2, [a, b])")
        # 2 is too close to 1 with a,b
        self.expect(results, ["GENE_3"])

        results = self.run_test("A", 10, 20, "c and not minimum(1, [a, b])")
        expected = {}  # 3 has an a, so it's out with a min of 1
        assert results == expected

        results = self.run_test("A", 50, 50, "c and not minimum(2, [a, b])")
        # even with min of 2, 3 now reaches 1, so there's another b
        assert results == {}

    def test_negated_group(self):
        # only gene 3 is far enough from both an a and an f
        results = self.run_test("A", 10, 10, "c and not (a or f)")
        self.expect(results, ["GENE_3"])

    def test_negated_cds_a(self):
        results = self.run_test("A", 10, 20, "c and not cds(a and b)")
        assert "GENE_1" not in results  # since 1 has both a and b
        assert "GENE_2" not in results  # since 1 is within range and has both

    def test_negated_cds_b(self):
        results = self.run_test("A", 5, 5, "c and not cds(a and c)")
        assert "GENE_2" not in results  # since 2 has a and c
        assert "GENE_3" in results  # since 3 has c but isn't in range of 2

    def test_negated_cds_c(self):
        results = self.run_test("A", 10, 20, "e and not cds(a or f)")
        # everything is close to a CDS with a or f, so no hits
        self.expect(results, [])

    def test_negated_cds_d(self):
        # GENE_1 excludes itself and everything in range
        # GENE_4+ excluded because they have no c in range
        results = self.run_test("A", 15, 20, "c and not cds(a and b)")
        self.expect(results, ["GENE_3"])

    def test_cds_single(self):
        with self.assertRaisesRegex(rule_parser.RuleSyntaxError, "must contain more than a single identifier"):
            self.run_test("A", 15, 20, "cds(a)")

    def test_negated_singles(self):
        results = self.run_test("A", 10, 20, "c and not a and not e and not f")
        self.expect(results, ["GENE_3"])

        results = self.run_test("A", 10, 20, "f and not a and not b and not c")
        self.expect(results, ["GENE_4", "GENE_5"])

    def test_minscore(self):
        # higher than default, will not recognise GENE_4
        results = self.run_test("A", 10, 20, "minscore(e, 100)")
        self.expect(results, [])

        # on the cutoff (> vs >=)
        results = self.run_test("A", 10, 20, "minscore(e, 50)")
        self.expect(results, ["GENE_4"])

        # below the cutoff
        results = self.run_test("A", 10, 20, "minscore(e, 49)")
        print(self.results_by_id["GENE_4"])
        self.expect(results, ["GENE_4"])

    def test_negated_minscore(self):
        # higher than default, everything with f counts
        results = self.run_test("A", 1, 1, "f and not minscore(e, 100)")
        self.expect(results, ["GENE_4", "GENE_5"])

        # on the cutoff (> vs >=), still no hit
        results = self.run_test("A", 1, 1, "f and not minscore(e, 50)")
        self.expect(results, ["GENE_5"])

        # below the cutoff
        results = self.run_test("A", 1, 1, "f and not minscore(e, 49)")
        self.expect(results, ["GENE_5"])


class RuleParserTest(unittest.TestCase):
    def setUp(self):
        self.signature_names = set(["a", "b", "c", "d", "more", "other"])

    def parse(self, text):
        return rule_parser.Parser(text, self.signature_names)

    def test_invalid_signature(self):
        with self.assertRaises(ValueError) as details:
            self.parse(format_as_rule("A", 10, 20, "badname or a"))
        assert str(details.exception) == "Rules contained identifers without signatures: badname"

    def test_stringify(self):
        rule_lines = [format_as_rule(*args) for args in [["abc", 10, 20, "a and b or not (c and d)"],
                                                         ["def", 7, 30, "minimum(3, [a, b, other]) and more"],
                                                         ["fgh", 15, 20, "c and not cds(a and b)"]]]
        rules = self.parse("\n".join(rule_lines)).rules
        assert rule_lines == [rule.reconstruct_rule_text() for rule in rules]

    def test_extra_whitespace(self):
        rules = self.parse("RULE A     CUTOFF 10\tEXTENT\t20\n CONDITIONS \t  a").rules
        assert len(rules) == 1
        assert str(rules[0]) == "A\t10\t20\ta"

    def test_cutoff_extent_parsing(self):
        rules = self.parse(format_as_rule("A", 10, 20, "a")).rules
        assert len(rules) == 1
        assert rules[0].cutoff == 10000
        assert rules[0].extent == 20000

        with self.assertRaises(rule_parser.RuleSyntaxError):
            self.parse("A 10 a or b")
        with self.assertRaises(rule_parser.RuleSyntaxError):
            self.parse("A b 10 a or b")

        with self.assertRaises(rule_parser.RuleSyntaxError):
            self.parse("RULE A CUTOFF 10 CONDITIONS a or b")
        with self.assertRaises(rule_parser.RuleSyntaxError):
            self.parse("RULE A CUTOFF b EXTENT 10 CONDITIONS a or b")

    def test_inline_comments(self):
        rule_chunk = "RULE name CUTOFF 20 EXTENT 20 CONDITIONS a"
        rules = self.parse("# comment line\n" + rule_chunk).rules
        assert len(rules) == 1 and rules[0].name == "name"

        rules = self.parse(rule_chunk + "#comment").rules
        assert len(rules) == 1 and rules[0].name == "name"

    def test_section_comments(self):
        rules = self.parse("RULE name COMMENT this is a section comment CUTOFF 20"
                           " EXTENT 20 CONDITIONS a").rules
        assert len(rules) == 1 and rules[0].comments == "this is a section comment"

    def test_single_superior(self):
        rules = self.parse("RULE first CUTOFF 20 EXTENT 20 CONDITIONS a "
                           "RULE sub SUPERIORS first CUTOFF 20 EXTENT 20 CONDITIONS c"
                           ).rules_by_name
        assert len(rules) == 2
        assert rules["sub"].superiors == ["first"]
        assert rules["first"].superiors == []

    def test_multiple_superiors(self):
        rules = self.parse("RULE first CUTOFF 20 EXTENT 20 CONDITIONS a "
                           "RULE second CUTOFF 20 EXTENT 20 CONDITIONS b "
                           "RULE sub SUPERIORS first,second CUTOFF 20 EXTENT 20 CONDITIONS c"
                           ).rules_by_name
        assert len(rules) == 3 and rules["sub"].superiors == ["first", "second"]
        assert not rules["first"].superiors and not rules["second"].superiors

    def test_unknown_superiors(self):
        # bad ordering
        with self.assertRaisesRegex(ValueError, "Unknown rule name: second"):
            self.parse("RULE first CUTOFF 20 EXTENT 20 CONDITIONS a "
                       "RULE sub SUPERIORS first,second CUTOFF 20 EXTENT 20 CONDITIONS c"
                       "RULE second CUTOFF 20 EXTENT 20 CONDITIONS b ")
        # completely undefined
        with self.assertRaisesRegex(ValueError, "Unknown rule name: second"):
            self.parse("RULE first CUTOFF 20 EXTENT 20 CONDITIONS a "
                       "RULE sub SUPERIORS first,second CUTOFF 20 EXTENT 20 CONDITIONS c")

    def test_self_referential_superiors(self):
        with self.assertRaisesRegex(ValueError, "Unknown rule name: sub"):
            self.parse("RULE first CUTOFF 20 EXTENT 20 CONDITIONS a "
                       "RULE sub SUPERIORS sub,first CUTOFF 20 EXTENT 20 CONDITIONS c")

    def test_empty_superiors(self):
        with self.assertRaisesRegex(rule_parser.RuleSyntaxError,
                                    "Expected identifier but found cutoff"):
            self.parse("RULE first CUTOFF 20 EXTENT 20 CONDITIONS a "
                       "RULE sub SUPERIORS CUTOFF 20 EXTENT 20 CONDITIONS c")

    def test_chained_superiors(self):
        with self.assertRaisesRegex(ValueError, "A rule cannot have a superior which has its own superior"):
            self.parse("RULE first CUTOFF 20 EXTENT 20 CONDITIONS a "
                       "RULE second SUPERIORS first CUTOFF 20 EXTENT 20 CONDITIONS b "
                       "RULE sub SUPERIORS second CUTOFF 20 EXTENT 20 CONDITIONS c")

    def test_related(self):
        rules = self.parse("RULE name RELATED b, c CUTOFF 20 EXTENT 20 CONDITIONS a").rules
        assert rules[0].related == ["b", "c"]

    def test_empty_related(self):
        with self.assertRaises(rule_parser.RuleSyntaxError):
            self.parse("RULE name RELATED CUTOFF 20 EXTENT 20 CONDITIONS a")

        with self.assertRaises(rule_parser.RuleSyntaxError):
            self.parse("RULE name RELATED")

    def test_missing_group_close(self):
        with self.assertRaises(rule_parser.RuleSyntaxError):
            self.parse(format_as_rule("A", 10, 10, "(a or b"))

    def test_missing_group_open(self):
        with self.assertRaises(rule_parser.RuleSyntaxError):
            self.parse(format_as_rule("A", 10, 10, "a or b) and c"))

    def test_bad_minscore(self):
        # negative scores
        with self.assertRaises(ValueError):
            self.parse(format_as_rule("A", 10, 20, "minscore(e, -1)"))

        # missing/invalid score syntax
        with self.assertRaises(rule_parser.RuleSyntaxError):
            self.parse(format_as_rule("A", 10, 20, "minscore(e)"))
        with self.assertRaises(rule_parser.RuleSyntaxError):
            self.parse(format_as_rule("A", 10, 20, "minscore(e,)"))
        with self.assertRaises(rule_parser.RuleSyntaxError):
            self.parse(format_as_rule("A", 10, 20, "minscore(e -1)"))

        # missing identifer
        with self.assertRaises(rule_parser.RuleSyntaxError):
            self.parse(format_as_rule("A", 10, 20, "minscore(5)"))
        with self.assertRaises(rule_parser.RuleSyntaxError):
            self.parse(format_as_rule("A", 10, 20, "minscore(,5)"))

        # empty
        with self.assertRaises(rule_parser.RuleSyntaxError):
            self.parse(format_as_rule("A", 10, 20, "minscore()"))

    def test_missing_binary_operand(self):
        with self.assertRaises(rule_parser.RuleSyntaxError):
            self.parse(format_as_rule("A", 10, 10, "(a or )"))

    def test_missing_op(self):
        with self.assertRaises(rule_parser.RuleSyntaxError):
            self.parse(format_as_rule("A", 10, 10, "(a  b)"))

    def test_double_not(self):
        with self.assertRaises(rule_parser.RuleSyntaxError):
            self.parse(format_as_rule("A", 10, 10, "not not a"))

    def test_only_not(self):
        with self.assertRaises(rule_parser.RuleSyntaxError):
            self.parse(format_as_rule("A", 10, 10, "not"))

    def test_repeated_minimum(self):
        with self.assertRaises(ValueError):
            self.parse(format_as_rule("A", 10, 10, "minimum(2, [a, a])"))

    def test_empty_minimum(self):
        with self.assertRaises(rule_parser.RuleSyntaxError):
            self.parse(format_as_rule("A", 10, 10, "minimum()"))

    def test_minimum_without_group(self):
        with self.assertRaises(rule_parser.RuleSyntaxError):
            self.parse(format_as_rule("A", 10, 10, "minimum a"))

    def test_minimum_with_bad_count(self):
        with self.assertRaises(rule_parser.RuleSyntaxError):
            self.parse(format_as_rule("A", 10, 10, "minimum([a,b])"))
        with self.assertRaises(rule_parser.RuleSyntaxError):
            self.parse(format_as_rule("A", 10, 10, "minimum(1.3, [a,b])"))
        with self.assertRaises(rule_parser.RuleSyntaxError):
            self.parse(format_as_rule("A", 10, 10, "minimum(1.3 [a,b])"))
        with self.assertRaises(ValueError):
            self.parse(format_as_rule("A", 10, 10, "minimum(0, [a,b])"))
        with self.assertRaises(ValueError):
            self.parse(format_as_rule("A", 10, 10, "minimum(-3, [a,b])"))

    def test_bad_identifier(self):
        assert not rule_parser.is_legal_identifier("a.b")
        with self.assertRaises(rule_parser.RuleSyntaxError):
            self.parse(format_as_rule("A", 10, 10, "a.b or c"))

        assert not rule_parser.is_legal_identifier("0sdf")
        with self.assertRaises(rule_parser.RuleSyntaxError):
            self.parse(format_as_rule("A", 10, 10, "a or 0sdf"))

        with self.assertRaises(rule_parser.RuleSyntaxError):
            self.parse(format_as_rule("A", 10, 10, "0sdf or a"))

        assert not rule_parser.is_legal_identifier("a!b")
        with self.assertRaises(rule_parser.RuleSyntaxError):
            self.parse(format_as_rule("A", 10, 10, "a or a!b"))

    def test_bad_syntax(self):
        with self.assertRaises(rule_parser.RuleSyntaxError):
            self.parse(format_as_rule("A", 10, 10, "name or name2(missing and op)"))
        with self.assertRaises(rule_parser.RuleSyntaxError):
            self.parse(format_as_rule("A", 10, 10, "name not name2"))
        with self.assertRaises(rule_parser.RuleSyntaxError):
            self.parse(format_as_rule("A", 10, 10, "name not or name2"))

    def test_bad_cds(self):
        with self.assertRaises(rule_parser.RuleSyntaxError):
            self.parse(format_as_rule("A", 10, 10, "cds()"))

        with self.assertRaises(rule_parser.RuleSyntaxError):
            self.parse(format_as_rule("A", 10, 10, "cds(1)"))

    def test_repetitions(self):
        with self.assertRaises(ValueError):
            self.parse(format_as_rule("A", 10, 10, "a or b or a"))
        with self.assertRaises(ValueError):
            self.parse(format_as_rule("A", 10, 10, "cds(a and b) or cds(a and b)"))

        # would be nice to catch these as redundant, but in the meantime
        self.parse(format_as_rule("A", 10, 10, "cds(a and b) or a and b"))
        self.parse(format_as_rule("A", 10, 10, "a and b or a"))

    def test_emptylines(self):
        rules = self.parse("\nRULE name CUTOFF 20 EXTENT 20 CONDITIONS a").rules
        assert len(rules) == 1 and rules[0].name == "name"

    def test_single_no_positive(self):
        rules = format_as_rule("A", 10, 10, "not a")
        with self.assertRaisesRegex(ValueError, "at least one positive requirement"):
            self.parse(rules)

    def test_and_no_positive(self):
        rules = format_as_rule("A", 10, 10, "not a and not c")
        with self.assertRaisesRegex(ValueError, "at least one positive requirement"):
            self.parse(rules)

    def test_or_no_positive(self):
        rules = format_as_rule("A", 10, 10, "not a or not c")
        with self.assertRaisesRegex(ValueError, "at least one positive requirement"):
            self.parse(rules)

    def test_cds_no_positive(self):
        rules = format_as_rule("A", 10, 10, "not cds(a and b)")
        with self.assertRaisesRegex(ValueError, "at least one positive requirement"):
            self.parse(rules)

        rules = format_as_rule("A", 10, 10, "not cds(a or b)")
        with self.assertRaisesRegex(ValueError, "at least one positive requirement"):
            self.parse(rules)

    def test_minimum_no_positive(self):
        rules = format_as_rule("A", 10, 10, "not minimum(2, [a, b])")
        with self.assertRaisesRegex(ValueError, "at least one positive requirement"):
            self.parse(rules)

    def test_score_no_positive(self):
        rules = format_as_rule("A", 10, 10, "not minscore(a, 15)")
        with self.assertRaisesRegex(ValueError, "at least one positive requirement"):
            self.parse(rules)

    def test_deep_no_positive(self):
        for rule in ["(not a) and (not b)", "not (a and b)",
                     "not a or (not b and not (c or d))"]:
            with self.assertRaisesRegex(ValueError, "at least one positive requirement"):
                self.parse(format_as_rule("A", 10, 10, rule))


class TokenTest(unittest.TestCase):
    def test_alpha(self):
        token = rule_parser.Token("abc", 1, 7)
        assert token.identifier == "abc"
        assert token.position == 7 - len(token.identifier)
        assert repr(token) == "'abc'"
        with self.assertRaises(AttributeError):
            # if it isn't an error, print it out for debugging
            print(token.value)

    def test_numeric(self):
        token = rule_parser.Token("12", 1, 6)
        assert token.value == 12
        assert token.position == 6 - len("12")
        assert repr(token) == ("12")
        with self.assertRaises(AttributeError):
            # if it isn't an error, print it out for debugging
            print(token.identifier)

    def test_operators(self):
        token = rule_parser.Token("minimum", 1, 20)
        assert token.type == rule_parser.TokenTypes.MINIMUM
        token = rule_parser.Token("and", 1, 20)
        assert token.type == rule_parser.TokenTypes.AND
        token = rule_parser.Token("or", 1, 20)
        assert token.type == rule_parser.TokenTypes.OR
        token = rule_parser.Token("not", 1, 20)
        assert token.type == rule_parser.TokenTypes.NOT

        # as a bonus, test the type converts to strings nicely for debugging
        assert str(token.type) == "not"
