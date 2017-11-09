# License: GNU Affero General Public License v3 or later
# A copy of GNU AGPL v3 should have been included in this software package in LICENSE.txt.

# for test files, silence irrelevant and noisy pylint warnings
# pylint: disable=no-self-use,protected-access,missing-docstring,too-many-public-methods

import unittest
from minimock import mock, restore

from antismash.common.test.helpers import DummyCDS, DummyRecord
from antismash.detection.hmm_detection import rule_parser
from antismash.detection.hmm_detection import signatures  # for mocking # pylint: disable=unused-import
from antismash.detection.hmm_detection.test.test_hmm_detection import FakeHSP


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
            "GENE_1": [FakeHSP("a", "GENE_1", 0, 10, 50, 0),
                       FakeHSP("b", "GENE_1", 0, 10, 50, 0)],
            "GENE_2": [FakeHSP("a", "GENE_1", 0, 10, 50, 0),
                       FakeHSP("c", "GENE_1", 0, 10, 50, 0)],
            "GENE_3": [FakeHSP("b", "GENE_1", 0, 10, 50, 0),
                       FakeHSP("c", "GENE_1", 0, 10, 50, 0)],
            "GENE_4": [FakeHSP("e", "GENE_1", 0, 10, 50, 0),
                       FakeHSP("f", "GENE_1", 0, 10, 50, 0)],
            "GENE_5": [FakeHSP("f", "GENE_1", 0, 10, 50, 0),
                       FakeHSP("g", "GENE_1", 0, 10, 50, 0)]}
        test_names = set(["a", "b", "c", "d", "e", "f", "g", "modelA", "modelB"])
        mock('signatures.get_signature_names', returns=test_names)

    def tearDown(self):
        restore()

    def run_test(self, name, cutoff, extent, conditions):
        rules = format_as_rule(name, cutoff, extent, conditions)
        rules = rule_parser.Parser(rules).rules
        for rule in rules:
            assert rule.contains_positive_condition()

        detected_types = {}
        cds_with_hits = sorted(self.results_by_id, key=lambda gene_id: self.feature_by_id[gene_id].location.start)
        for cds in cds_with_hits:
            detected_type = detected_types.get(cds)
            if detected_type:
                continue
            rule_results = [rule.detect(cds, self.feature_by_id, self.results_by_id) for rule in rules]
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
                FakeHSP("modelA", "GENE_1", 0, 10, 50, 0),
                FakeHSP("modelB", "GENE_1", 0, 10, 50, 0)
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
        test_names = set(["a", "b", "c", "d", "more", "other"])
        mock('signatures.get_signature_names', returns=test_names)

    def tearDown(self):
        restore()

    def test_invalid_signature(self):
        with self.assertRaises(ValueError) as details:
            rule_parser.Parser(format_as_rule("A", 10, 20, "badname or a"))
        assert str(details.exception) == "Rules contained identifers without signatures: badname"

    def test_stringify(self):
        rule_lines = [format_as_rule(*args) for args in [["abc", 10, 20, "a and b or not (c and d)"],
                                                         ["def", 7, 30, "minimum(3, [a, b, other]) and more"],
                                                         ["fgh", 15, 20, "c and not cds(a and b)"]]]
        rules = rule_parser.Parser("\n".join(rule_lines)).rules
        assert rule_lines == [rule.reconstruct_rule_text() for rule in rules]

    def test_extra_whitespace(self):
        rules = rule_parser.Parser("RULE A     CUTOFF 10\tEXTENT\t20\n CONDITIONS \t  a").rules
        assert len(rules) == 1
        assert str(rules[0]) == "A\t10\t20\ta"

    def test_cutoff_extent_parsing(self):
        rules = rule_parser.Parser(format_as_rule("A", 10, 20, "a")).rules
        assert len(rules) == 1
        assert rules[0].cutoff == 10000
        assert rules[0].extent == 20000

        with self.assertRaises(rule_parser.RuleSyntaxError):
            rule_parser.Parser("A 10 a or b")
        with self.assertRaises(rule_parser.RuleSyntaxError):
            rule_parser.Parser("A b 10 a or b")

        with self.assertRaises(rule_parser.RuleSyntaxError):
            rule_parser.Parser("RULE A CUTOFF 10 CONDITIONS a or b")
        with self.assertRaises(rule_parser.RuleSyntaxError):
            rule_parser.Parser("RULE A CUTOFF b EXTENT 10 CONDITIONS a or b")

    def test_comments(self):
        rule_chunk = "RULE name CUTOFF 20 EXTENT 20 CONDITIONS a"
        rules = rule_parser.Parser("# comment line\n" + rule_chunk).rules
        assert len(rules) == 1 and rules[0].name == "name"

        rules = rule_parser.Parser(rule_chunk + "#comment").rules
        assert len(rules) == 1 and rules[0].name == "name"

    def test_missing_group_close(self):
        with self.assertRaises(rule_parser.RuleSyntaxError):
            rule_parser.Parser(format_as_rule("A", 10, 10, "(a or b"))

    def test_missing_group_open(self):
        with self.assertRaises(rule_parser.RuleSyntaxError):
            rule_parser.Parser(format_as_rule("A", 10, 10, "a or b) and c"))

    def test_bad_minscore(self):
        # negative scores
        with self.assertRaises(ValueError):
            rule_parser.Parser(format_as_rule("A", 10, 20, "minscore(e, -1)"))

        # missing/invalid score syntax
        with self.assertRaises(rule_parser.RuleSyntaxError):
            rule_parser.Parser(format_as_rule("A", 10, 20, "minscore(e)"))
        with self.assertRaises(rule_parser.RuleSyntaxError):
            rule_parser.Parser(format_as_rule("A", 10, 20, "minscore(e,)"))
        with self.assertRaises(rule_parser.RuleSyntaxError):
            rule_parser.Parser(format_as_rule("A", 10, 20, "minscore(e -1)"))

        # missing identifer
        with self.assertRaises(rule_parser.RuleSyntaxError):
            rule_parser.Parser(format_as_rule("A", 10, 20, "minscore(5)"))
        with self.assertRaises(rule_parser.RuleSyntaxError):
            rule_parser.Parser(format_as_rule("A", 10, 20, "minscore(,5)"))

        # empty
        with self.assertRaises(rule_parser.RuleSyntaxError):
            rule_parser.Parser(format_as_rule("A", 10, 20, "minscore()"))

    def test_missing_binary_operand(self):
        with self.assertRaises(rule_parser.RuleSyntaxError):
            rule_parser.Parser(format_as_rule("A", 10, 10, "(a or )"))

    def test_missing_op(self):
        with self.assertRaises(rule_parser.RuleSyntaxError):
            rule_parser.Parser(format_as_rule("A", 10, 10, "(a  b)"))

    def test_double_not(self):
        with self.assertRaises(rule_parser.RuleSyntaxError):
            rule_parser.Parser(format_as_rule("A", 10, 10, "not not a"))

    def test_only_not(self):
        with self.assertRaises(rule_parser.RuleSyntaxError):
            rule_parser.Parser(format_as_rule("A", 10, 10, "not"))

    def test_repeated_minimum(self):
        with self.assertRaises(ValueError):
            rule_parser.Parser(format_as_rule("A", 10, 10, "minimum(2, [a, a])"))

    def test_empty_minimum(self):
        with self.assertRaises(rule_parser.RuleSyntaxError):
            rule_parser.Parser(format_as_rule("A", 10, 10, "minimum()"))

    def test_minimum_without_group(self):
        with self.assertRaises(rule_parser.RuleSyntaxError):
            rule_parser.Parser(format_as_rule("A", 10, 10, "minimum a"))

    def test_minimum_with_bad_count(self):
        with self.assertRaises(rule_parser.RuleSyntaxError):
            rule_parser.Parser(format_as_rule("A", 10, 10, "minimum([a,b])"))
        with self.assertRaises(rule_parser.RuleSyntaxError):
            rule_parser.Parser(format_as_rule("A", 10, 10, "minimum(1.3, [a,b])"))
        with self.assertRaises(rule_parser.RuleSyntaxError):
            rule_parser.Parser(format_as_rule("A", 10, 10, "minimum(1.3 [a,b])"))
        with self.assertRaises(ValueError):
            rule_parser.Parser(format_as_rule("A", 10, 10, "minimum(0, [a,b])"))
        with self.assertRaises(ValueError):
            rule_parser.Parser(format_as_rule("A", 10, 10, "minimum(-3, [a,b])"))

    def test_bad_identifier(self):
        assert not rule_parser.is_legal_identifier("a.b")
        with self.assertRaises(rule_parser.RuleSyntaxError):
            rule_parser.Parser(format_as_rule("A", 10, 10, "a.b or c"))

        assert not rule_parser.is_legal_identifier("0sdf")
        with self.assertRaises(rule_parser.RuleSyntaxError):
            rule_parser.Parser(format_as_rule("A", 10, 10, "a or 0sdf"))

        with self.assertRaises(rule_parser.RuleSyntaxError):
            rule_parser.Parser(format_as_rule("A", 10, 10, "0sdf or a"))

        assert not rule_parser.is_legal_identifier("a!b")
        with self.assertRaises(rule_parser.RuleSyntaxError):
            rule_parser.Parser(format_as_rule("A", 10, 10, "a or a!b"))

    def test_bad_syntax(self):
        with self.assertRaises(rule_parser.RuleSyntaxError):
            rule_parser.Parser(format_as_rule("A", 10, 10, "name or name2(missing and op)"))
        with self.assertRaises(rule_parser.RuleSyntaxError):
            rule_parser.Parser(format_as_rule("A", 10, 10, "name not name2"))
        with self.assertRaises(rule_parser.RuleSyntaxError):
            rule_parser.Parser(format_as_rule("A", 10, 10, "name not or name2"))

    def test_bad_cds(self):
        with self.assertRaises(rule_parser.RuleSyntaxError):
            rule_parser.Parser(format_as_rule("A", 10, 10, "cds()"))

        with self.assertRaises(rule_parser.RuleSyntaxError):
            rule_parser.Parser(format_as_rule("A", 10, 10, "cds(1)"))

    def test_repetitions(self):
        with self.assertRaises(ValueError):
            rule_parser.Parser(format_as_rule("A", 10, 10, "a or b or a"))
        with self.assertRaises(ValueError):
            rule_parser.Parser(format_as_rule("A", 10, 10, "cds(a and b) or cds(a and b)"))

        # would be nice to catch these as redundant, but in the meantime
        rule_parser.Parser(format_as_rule("A", 10, 10, "cds(a and b) or a and b"))
        rule_parser.Parser(format_as_rule("A", 10, 10, "a and b or a"))

    def test_emptylines(self):
        rules = rule_parser.Parser("\nRULE name CUTOFF 20 EXTENT 20 CONDITIONS a").rules
        assert len(rules) == 1 and rules[0].name == "name"

    def test_single_no_positive(self):
        rules = format_as_rule("A", 10, 10, "not a")
        with self.assertRaisesRegex(ValueError, "at least one positive requirement"):
            rule_parser.Parser(rules)

    def test_and_no_positive(self):
        rules = format_as_rule("A", 10, 10, "not a and not c")
        with self.assertRaisesRegex(ValueError, "at least one positive requirement"):
            rule_parser.Parser(rules)

    def test_or_no_positive(self):
        rules = format_as_rule("A", 10, 10, "not a or not c")
        with self.assertRaisesRegex(ValueError, "at least one positive requirement"):
            rule_parser.Parser(rules)

    def test_cds_no_positive(self):
        rules = format_as_rule("A", 10, 10, "not cds(a and b)")
        with self.assertRaisesRegex(ValueError, "at least one positive requirement"):
            rule_parser.Parser(rules)

        rules = format_as_rule("A", 10, 10, "not cds(a or b)")
        with self.assertRaisesRegex(ValueError, "at least one positive requirement"):
            rule_parser.Parser(rules)

    def test_minimum_no_positive(self):
        rules = format_as_rule("A", 10, 10, "not minimum(2, [a, b])")
        with self.assertRaisesRegex(ValueError, "at least one positive requirement"):
            rule_parser.Parser(rules)

    def test_score_no_positive(self):
        rules = format_as_rule("A", 10, 10, "not minscore(a, 15)")
        with self.assertRaisesRegex(ValueError, "at least one positive requirement"):
            rule_parser.Parser(rules)

    def test_deep_no_positive(self):
        for rule in ["(not a) and (not b)", "not (a and b)",
                     "not a or (not b and not (c or d))"]:
            with self.assertRaisesRegex(ValueError, "at least one positive requirement"):
                rule_parser.Parser(format_as_rule("A", 10, 10, rule))


class TokenTest(unittest.TestCase):
    def test_alpha(self):
        t = rule_parser.Token("abc", 1, 7)
        assert t.identifier == "abc"
        assert t.position == 7 - len(t.identifier)
        assert repr(t) == "'abc'"
        with self.assertRaises(AttributeError):
            # if it isn't an error, print it out for debugging
            print(t.value)

    def test_numeric(self):
        t = rule_parser.Token("12", 1, 6)
        assert t.value == 12
        assert t.position == 6 - len("12")
        assert repr(t) == ("12")
        with self.assertRaises(AttributeError):
            # if it isn't an error, print it out for debugging
            print(t.identifier)

    def test_operators(self):
        t = rule_parser.Token("minimum", 1, 20)
        assert t.type == rule_parser.TokenTypes.MINIMUM
        t = rule_parser.Token("and", 1, 20)
        assert t.type == rule_parser.TokenTypes.AND
        t = rule_parser.Token("or", 1, 20)
        assert t.type == rule_parser.TokenTypes.OR
        t = rule_parser.Token("not", 1, 20)
        assert t.type == rule_parser.TokenTypes.NOT

        # as a bonus, test the type converts to strings nicely for debugging
        assert str(t.type) == "not"
