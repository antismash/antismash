# License: GNU Affero General Public License v3 or later
# A copy of GNU AGPL v3 should have been included in this software package in LICENSE.txt.

import unittest
from minimock import mock, restore

import antismash.modules.hmm_detection.rule_parser as rule_parser
from antismash.modules.hmm_detection import signatures  # for mocking
from antismash.modules.hmm_detection.test.test_hmm_detection import FakeFeature, FakeHSP, FakeRecord, FeatureLocation

class DetectionTest(unittest.TestCase):
    def setUp(self):
        self.feature_by_id = {
            "GENE_1" : FakeFeature("CDS", FeatureLocation(0, 30000), {"locus_tag": ["GENE_1"]}),
            "GENE_2" : FakeFeature("CDS", FeatureLocation(30000, 50000), {"locus_tag": ["GENE_2"]}),
            "GENE_3" : FakeFeature("CDS", FeatureLocation(70000, 90000), {"locus_tag": ["GENE_3"]}),
            "GENE_X" : FakeFeature("CDS", FeatureLocation(95000, 100000), {"locus_tag": ["GENE_X"]}),
            "GENE_4" : FakeFeature("CDS", FeatureLocation(125000, 140000), {"locus_tag": ["GENE_4"]}),
            "GENE_5" : FakeFeature("CDS", FeatureLocation(145000, 150000), {"locus_tag": ["GENE_5"]})
        }
        self.features = []
        for gene_id in self.feature_by_id:
            self.features.append(self.feature_by_id[gene_id])
        self.features.sort(key=lambda x: x.location.start) # vital for py3 < 3.5
        self.record = FakeRecord(self.features)

        self.results_by_id = {
            "GENE_1" : [FakeHSP("a", "GENE_1", 0, 10, 50, 0),
                        FakeHSP("b", "GENE_1", 0, 10, 50, 0)],
            "GENE_2" : [FakeHSP("a", "GENE_1", 0, 10, 50, 0),
                        FakeHSP("c", "GENE_1", 0, 10, 50, 0)],
            "GENE_3" : [FakeHSP("b", "GENE_1", 0, 10, 50, 0),
                        FakeHSP("c", "GENE_1", 0, 10, 50, 0)],
            "GENE_4" : [FakeHSP("e", "GENE_1", 0, 10, 50, 0),
                        FakeHSP("f", "GENE_1", 0, 10, 50, 0)],
            "GENE_5" : [FakeHSP("f", "GENE_1", 0, 10, 50, 0),
                        FakeHSP("g", "GENE_1", 0, 10, 50, 0)]}
        test_names = set(["a", "b", "c", "d", "e", "f", "g", "modelA", "modelB"])
        mock('signatures.get_signature_names', returns=test_names)

    def tearDown(self):
        restore()

    def run_test(self, rules):
        enabled_clustertypes = [rule.split()[0] for rule in rules]
        rules = rule_parser.Parser(rules).rules
        rules = [rule for rule in rules if rule.name in enabled_clustertypes]

        detected_types = {}
        cds_with_hits = sorted(self.results_by_id, key = lambda gene_id: self.feature_by_id[gene_id].location.start)
        for cds in cds_with_hits:
            detected_type = detected_types.get(cds)
            if detected_type:
                continue
            results = [rule.name for rule in rules if rule.detect(cds, self.feature_by_id, self.results_by_id)]
            if results:
                detected_types[cds] = results

        for gid in detected_types:
            detected_types[gid] = set(detected_types[gid])
        return detected_types

    def expect(self, results, genes_to_hit):
        expected = {gene : set("A") for gene in genes_to_hit}
        assert results == expected

    def test_single(self):
        results = self.run_test(["A 10 20 a"])
        self.expect(results, ["GENE_1", "GENE_2"])

    def test_simple_or(self):
        results = self.run_test(["A 10 20 a or c"])
        self.expect(results, ["GENE_1", "GENE_2", "GENE_3"])

    def test_chained_or(self):
        results = self.run_test(["A 10 20 a or b or c or e or f"])
        assert len(results) == 5 # all valid

    def test_simple_and(self):
        results = self.run_test(["A 10 20 a and c"])
        self.expect(results, ["GENE_1", "GENE_2"])

        results = self.run_test(["A 10 20 e and c"])
        self.expect(results, [])

    def test_cds_and(self):
        results = self.run_test(["A 10 20 cds(a and c)"])
        # GENE_2 is an exact match
        # GENE_1 is in range, but doesn't satisfy any condition locally
        self.expect(results, ["GENE_2"])

    def test_cds_or(self):
        results = self.run_test(["A 10 20 cds(a or c)"])
        self.expect(results, ["GENE_1", "GENE_2", "GENE_3"])

    def test_cds_internal_combo(self):
        results = self.run_test(["A 10 20 cds(a and (b or c))"])
        self.expect(results, ["GENE_1", "GENE_2"])

    def test_chained_and(self):
        results = self.run_test(["A 10 20 a and b and not c"])
        self.expect(results, ["GENE_1"]) # 2 reaches c

        results = self.run_test(["A 10 20 a and b and not cds(a and b)"])
        # the range is too short to reach 2 for the a
        self.expect(results, [])

        results = self.run_test(["A 21 20 a and b and not cds(a and b)"])
        self.expect(results, ["GENE_3"]) # 3 has b, 2 has a, 1 reaches 2 with a,b

    def test_simple_minimum(self):
        results = self.run_test(["A 10 20 minimum(2, [a, b])"])
        self.expect(results, ["GENE_1", "GENE_2"])

        results = self.run_test(["A 100 20 minimum(2, [b, e])"])
        self.expect(results, ["GENE_1", "GENE_3", "GENE_4"])

    def test_single_gene(self):
        self.results_by_id = {
            "GENE_1" : [
                FakeHSP("modelA", "GENE_1", 0, 10, 50, 0),
                FakeHSP("modelB", "GENE_1", 0, 10, 50, 0)
            ]}
        self.feature_by_id = {
            "GENE_1" : FakeFeature("CDS", FeatureLocation(0, 30000), {"locus_tag": ["GENE_1"]})
        }

        results = self.run_test(["A 10 20 minimum(2, [modelA,modelB])"])
        self.expect(results, ["GENE_1"])

        results = self.run_test(["A 10 20 cds(modelA and modelB)"])
        self.expect(results, ["GENE_1"])

        results = self.run_test(["A 10 20 cds(modelA or modelB)"])
        self.expect(results, ["GENE_1"])

        results = self.run_test(["A 10 20 modelA and modelB"])
        self.expect(results, ["GENE_1"])

        results = self.run_test(["A 10 20 modelA or modelB"])
        self.expect(results, ["GENE_1"])

    def test_negated_minimum(self):
        results = self.run_test(["A 10 20 c and not minimum(2, [a, b])"])
        self.expect(results, ["GENE_3"])

        results = self.run_test(["A 10 20 c and not minimum(1, [a, b])"])
        expected = {} # 3 had a,c and 2 contributes a second option, b
        assert results == expected

        results = self.run_test(["A 50 50 c and not minimum(2, [a, b])"])
        # 3 now reaches a so it's excluded, but 4 reaches a c
        self.expect(results, ["GENE_4"])

    def test_negated_group(self):
        #only gene 3 is far enough from both an a and an f
        results = self.run_test(["A 10 10 c and not (a or f)"])
        self.expect(results, ["GENE_3"])

    def test_negated_cds(self):
        results = self.run_test(["A 10 20 not cds(a and b)"])
        assert "GENE_1" not in results # since 1 has both a and b
        assert "GENE_2" not in results # since 1 is within range and has both

        results = self.run_test(["A 10 20 a and not cds(a and b)"])
        assert "GENE_1" not in results # since 1 has both a and b
        assert "GENE_2" not in results # since 1 is within range and has both

        results = self.run_test(["A 10 20 not cds(a or f)"])
        # everything is close to a CDS with a or f, so no hits
        self.expect(results, [])

        # GENE_1 excludes itself and everything in range
        # GENE_4+ excluded because they have no c in range
        results = self.run_test(["A 15 20 c and not cds(a and b)"])
        self.expect(results, ["GENE_3"])

    def test_negated_singles(self):
        results = self.run_test(["A 10 20 not a and not e and not f"])
        self.expect(results, ["GENE_3"])

        results = self.run_test(["A 10 20 not a and not b and not c"])
        self.expect(results, ["GENE_4", "GENE_5"])

    def test_minscore(self):
        # higher than default, will not recognise GENE_4
        results = self.run_test(["A 10 20 minscore(e, 100)"])
        self.expect(results, [])

        # on the cutoff (> vs >=), still no hit
        results = self.run_test(["A 10 20 minscore(e, 50)"])
        self.expect(results, ["GENE_4"])

        # below the cutoff
        results = self.run_test(["A 10 20 minscore(e, 49)"])
        print(self.results_by_id["GENE_4"])
        self.expect(results, ["GENE_4"])


    def test_negated_minscore(self):
        # higher than default, everything with f counts
        results = self.run_test(["A 1 1 f and not minscore(e, 100)"])
        self.expect(results, ["GENE_4", "GENE_5"])

        # on the cutoff (> vs >=), still no hit
        results = self.run_test(["A 1 1 f and not minscore(e, 50)"])
        self.expect(results, ["GENE_5"])

        # below the cutoff
        results = self.run_test(["A 1 1 f and not minscore(e, 49)"])
        self.expect(results, ["GENE_5"])

class RuleParserTest(unittest.TestCase):
    def setUp(self):
        test_names = set(["a", "b", "c", "d", "more", "other"])
        mock('signatures.get_signature_names', returns=test_names)

    def tearDown(self):
        restore()

    def test_invalid_signature(self):
        with self.assertRaises(ValueError) as details:
            rule_parser.Parser(["A 10 20 badname or a"])
        assert str(details.exception) == "Rules contained identifers without signatures: badname"

    def test_stringify(self):
        rule_lines = ["abc\t10\t20\ta and b or not (c and d)",
                      "def\t7\t30\tminimum(3, [a, b, other]) and more",
                      "fgh\t15\t20\tc and not cds(a and b)"]
        rules = rule_parser.Parser(rule_lines).rules
        assert rule_lines == list(map(str, rules))

    def test_extra_whitespace(self):
        rules = rule_parser.Parser(["A      10\t20  \t  a"]).rules
        assert len(rules) == 1
        assert str(rules[0]) == "A\t10\t20\ta"

    def test_cutoff_extent_parsing(self):
        rules = rule_parser.Parser(["A 10 20 a"]).rules
        assert len(rules) == 1
        assert rules[0].cutoff == 10000
        assert rules[0].extent == 20000

        with self.assertRaises(rule_parser.RuleSyntaxError):
            rule_parser.Parser(["A 10 a or b"])
        with self.assertRaises(rule_parser.RuleSyntaxError):
            rule_parser.Parser(["A b 10 a or b"])

    def test_comments(self):
        rules = rule_parser.Parser(["# comment line", "name 20 20 a"]).rules
        assert len(rules) == 1 and rules[0].name == "name"

        rules = rule_parser.Parser(["name 20 20 a #comment"]).rules
        assert len(rules) == 1 and rules[0].name == "name"

    def test_missing_group_close(self):
        with self.assertRaises(rule_parser.RuleSyntaxError):
            rule_parser.Parser(["A 10 10 (a or b"])

    def test_missing_group_open(self):
        with self.assertRaises(rule_parser.RuleSyntaxError):
            rule_parser.Parser(["A 10 10 a or b) and c"])

    def test_bad_minscore(self):
        # negative scores
        with self.assertRaises(ValueError):
            rule_parser.Parser(["A 10 20 minscore(e, -1)"])

        # missing/invalid score syntax
        with self.assertRaises(rule_parser.RuleSyntaxError):
            rule_parser.Parser(["A 10 20 minscore(e)"])
        with self.assertRaises(rule_parser.RuleSyntaxError):
            rule_parser.Parser(["A 10 20 minscore(e,)"])
        with self.assertRaises(rule_parser.RuleSyntaxError):
            rule_parser.Parser(["A 10 20 minscore(e -1)"])

        # missing identifer
        with self.assertRaises(rule_parser.RuleSyntaxError):
            rule_parser.Parser(["A 10 20 minscore(5)"])
        with self.assertRaises(rule_parser.RuleSyntaxError):
            rule_parser.Parser(["A 10 20 minscore(,5)"])

        # empty
        with self.assertRaises(rule_parser.RuleSyntaxError):
            rule_parser.Parser(["A 10 20 minscore()"])

    def test_missing_binary_operand(self):
        with self.assertRaises(rule_parser.RuleSyntaxError):
            rule_parser.Parser(["A 10 10 (a or )"])

    def test_missing_op(self):
        with self.assertRaises(rule_parser.RuleSyntaxError):
            rule_parser.Parser(["A 10 10 (a  b)"])

    def test_double_not(self):
        with self.assertRaises(rule_parser.RuleSyntaxError):
            rule_parser.Parser(["A 10 10 not not a"])

    def test_only_not(self):
        with self.assertRaises(rule_parser.RuleSyntaxError):
            rule_parser.Parser(["A 10 10 not"])

    def test_repeated_minimum(self):
        with self.assertRaises(ValueError):
            rule_parser.Parser(["A 10 10 minimum(2, [a, a])"])

    def test_empty_minimum(self):
        with self.assertRaises(rule_parser.RuleSyntaxError):
            rule_parser.Parser(["A 10 10 minimum()"])

    def test_minimum_without_group(self):
        with self.assertRaises(rule_parser.RuleSyntaxError):
            rule_parser.Parser(["A 10 10 minimum a"])

    def test_minimum_with_bad_count(self):
        with self.assertRaises(rule_parser.RuleSyntaxError):
            rule_parser.Parser(["A 10 10 minimum([a,b])"])
        with self.assertRaises(rule_parser.RuleSyntaxError):
            rule_parser.Parser(["A 10 10 minimum(1.3, [a,b])"])
        with self.assertRaises(rule_parser.RuleSyntaxError):
            rule_parser.Parser(["A 10 10 minimum(1.3 [a,b])"])
        with self.assertRaises(ValueError):
            rule_parser.Parser(["A 10 10 minimum(0, [a,b])"])
        with self.assertRaises(ValueError):
            rule_parser.Parser(["A 10 10 minimum(-3, [a,b])"])

    def test_bad_identifier(self):
        assert not rule_parser.is_legal_identifier("a.b")
        with self.assertRaises(rule_parser.RuleSyntaxError):
            rule_parser.Parser(["A 10 10 a.b or c"])

        assert not rule_parser.is_legal_identifier("0sdf")
        with self.assertRaises(rule_parser.RuleSyntaxError):
            rule_parser.Parser(["A 10 10 a or 0sdf"])

        with self.assertRaises(rule_parser.RuleSyntaxError):
            rule_parser.Parser(["A 10 10 0sdf or a"])

        assert not rule_parser.is_legal_identifier("a!b")
        with self.assertRaises(rule_parser.RuleSyntaxError):
            rule_parser.Parser(["A 10 10 a or a!b"])

    def test_bad_syntax(self):
        with self.assertRaises(rule_parser.RuleSyntaxError):
            rule_parser.Parser(["A 10 10 name or name2(missing and op)"])
        with self.assertRaises(rule_parser.RuleSyntaxError):
            rule_parser.Parser(["A 10 10 name not name2"])
        with self.assertRaises(rule_parser.RuleSyntaxError):
            rule_parser.Parser(["A 10 10 name not or name2"])

    def test_bad_cds(self):
        with self.assertRaises(rule_parser.RuleSyntaxError):
            rule_parser.Parser(["A 10 10 cds()"])

        with self.assertRaises(rule_parser.RuleSyntaxError):
            rule_parser.Parser(["A 10 10 cds(1)"])

    def test_repetitions(self):
        with self.assertRaises(ValueError):
            rule_parser.Parser(["A 10 10 a or b or a"])
        with self.assertRaises(ValueError):
            rule_parser.Parser(["A 10 10 cds(a and b) or cds(a and b)"])

        # would be nice to catch these as redundant, but in the meantime
        rule_parser.Parser(["A 10 10 cds(a and b) or a and b"])
        rule_parser.Parser(["A 10 10 a and b or a"])

    def test_emptylines(self):
        rules = rule_parser.Parser(["", "name 20 20 a"]).rules
        assert len(rules) == 1 and rules[0].name == "name"

class TokenTest(unittest.TestCase):
    def test_alpha(self):
        t = rule_parser.Token("abc", 7)
        assert t.identifier == "abc"
        assert t.position == 7 - len(t.identifier)
        assert repr(t) == "'abc'"
        with self.assertRaises(AttributeError):
            # if it isn't an error, print it out for debugging
            print(t.value)

    def test_numeric(self):
        t = rule_parser.Token("12", 6)
        assert t.value == 12
        assert t.position == 6 - len("12")
        assert repr(t) == ("12")
        with self.assertRaises(AttributeError):
            # if it isn't an error, print it out for debugging
            print(t.identifier)

    def test_operators(self):
        t = rule_parser.Token("minimum", 20)
        assert t.type == rule_parser.TokenTypes.MINIMUM
        t = rule_parser.Token("and", 20)
        assert t.type == rule_parser.TokenTypes.AND
        t = rule_parser.Token("or", 20)
        assert t.type == rule_parser.TokenTypes.OR
        t = rule_parser.Token("not", 20)
        assert t.type == rule_parser.TokenTypes.NOT

        # as a bonus, test the type converts to strings nicely for debugging
        assert str(t.type) == "not"
