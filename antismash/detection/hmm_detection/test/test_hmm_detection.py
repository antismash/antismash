# License: GNU Affero General Public License v3 or later
# A copy of GNU AGPL v3 should have been included in this software package in LICENSE.txt.

# for test files, silence irrelevant and noisy pylint warnings
# pylint: disable=use-implicit-booleaness-not-comparison,protected-access,missing-docstring,consider-using-with

from argparse import Namespace
import glob
import json
import importlib
import os
import pkgutil
import unittest
from unittest.mock import patch

from Bio.Seq import Seq

from antismash.common import path
from antismash.common.hmm_rule_parser import rule_parser, cluster_prediction as hmm_detection  # TODO: redo tests
from antismash.common.hmm_rule_parser.test.helpers import check_hmm_signatures, create_ruleset
from antismash.common.secmet import Record
from antismash.common.test.helpers import DummyCDS, DummyRecord, FakeHSPHit
from antismash.config import build_config, destroy_config
import antismash.detection.hmm_detection as core
from antismash.detection.hmm_detection import DynamicProfile, signatures


class HmmDetectionTest(unittest.TestCase):
    def setUp(self):
        self.config = build_config([])
        self.rules_file = path.get_full_path(__file__, "..", "cluster_rules", "strict.txt")
        self.signature_file = path.get_full_path(__file__, "..", "data", "hmmdetails.txt")
        self.signature_names = {sig.name for sig in core.get_signature_profiles()}.union(core.DYNAMIC_PROFILES)
        self.valid_categories = {cat.name for cat in core.get_rule_categories()}
        self.filter_file = path.get_full_path(__file__, "..", "filterhmmdetails.txt")
        self.results_by_id = {
            "GENE_1": [
                FakeHSPHit("modelA", "GENE_1", 0, 10, 50, 0),
                FakeHSPHit("modelB", "GENE_1", 0, 10, 50, 0)
            ],
            "GENE_2": [
                FakeHSPHit("modelC", "GENE_2", 0, 10, 50, 0),
                FakeHSPHit("modelB", "GENE_2", 0, 10, 50, 0)
            ],
            "GENE_3": [
                FakeHSPHit("modelC", "GENE_3", 0, 10, 50, 0),
                FakeHSPHit("modelF", "GENE_3", 0, 10, 50, 0)
            ],
            "GENE_4": [
                FakeHSPHit("modelA", "GENE_4", 0, 10, 50, 0),
                FakeHSPHit("modelE", "GENE_4", 0, 10, 50, 0)
            ],
            "GENE_5": [
                FakeHSPHit("modelA", "GENE_5", 0, 10, 50, 0),
                FakeHSPHit("modelG", "GENE_5", 0, 10, 50, 0)
            ]
        }
        self.feature_by_id = {
            "GENE_1": DummyCDS(0, 30000, locus_tag="GENE_1"),
            "GENE_2": DummyCDS(30000, 50000, locus_tag="GENE_2"),
            "GENE_3": DummyCDS(70000, 90000, locus_tag="GENE_3"),
            "GENE_X": DummyCDS(95000, 100000, locus_tag="GENE_X"),  # no hits
            "GENE_4": DummyCDS(125000, 140000, locus_tag="GENE_4"),
            "GENE_5": DummyCDS(130000, 150000, locus_tag="GENE_5")
        }

        self.test_names = {"modelA", "modelB", "modelC", "modelF", "modelG",
                           "a", "b", "c", "d"}

        self.categories = {"Cat"}

        self.rules = rule_parser.Parser("\n".join([
            "RULE MetaboliteA CATEGORY Cat CUTOFF 10 NEIGHBOURHOOD 5 CONDITIONS modelA",
            "RULE MetaboliteB CATEGORY Cat CUTOFF 10 NEIGHBOURHOOD 5 CONDITIONS cds(modelA and modelB)",
            "RULE MetaboliteC CATEGORY Cat CUTOFF 10 NEIGHBOURHOOD 5 CONDITIONS (modelA and modelB)",
            "RULE MetaboliteD CATEGORY Cat CUTOFF 20 NEIGHBOURHOOD 5 CONDITIONS minimum(2,[modelC,modelB]) and modelA",
            "RULE Metabolite0 CATEGORY Cat CUTOFF 1 NEIGHBOURHOOD 3 CONDITIONS modelF",
            "RULE Metabolite1 CATEGORY Cat CUTOFF 1 NEIGHBOURHOOD 3 CONDITIONS modelG"]),
            self.test_names, self.categories).rules
        self.record = Record()
        self.record._record.seq = Seq("A"*150000)
        for feature in self.feature_by_id.values():
            self.record.add_cds_feature(feature)

    def tearDown(self):
        # clear out any leftover config adjustments
        destroy_config()

    def test_overlaps_but_not_contains(self):
        # should get gene2 and gene3
        rules = rule_parser.Parser("\n".join([
                "RULE Overlap CATEGORY Cat CUTOFF 25 NEIGHBOURHOOD 5 CONDITIONS modelB and modelF "
                "RULE OverlapImpossible CATEGORY Cat CUTOFF 25 NEIGHBOURHOOD 5 CONDITIONS modelA and modelF"]),
                self.test_names, self.categories).rules
        detected_types, cluster_type_hits = hmm_detection.apply_cluster_rules(self.record, self.results_by_id, rules)
        assert detected_types == {"GENE_2": {"Overlap": {"modelB"}},
                                  "GENE_3": {"Overlap": {"modelF"}}}

        assert cluster_type_hits == {"Overlap": {"GENE_2", "GENE_3"}}

        # only 1 cluster should be found, since it requires both genes
        # if forming clusters by .is_contained_by(), 2 clusters will be formed
        # if finding rule hits uses .is_contained_by(), no clusters will be formed
        rules_by_name = {rule.name: rule for rule in rules}
        dummy_extender_domains = {}
        clusters = hmm_detection.find_protoclusters(self.record, cluster_type_hits, rules_by_name,
                                                    self.results_by_id, dummy_extender_domains)
        assert len(clusters) == 1
        assert clusters[0].product == "Overlap"
        assert clusters[0].core_location.start == 30000
        assert clusters[0].core_location.end == 90000

    def test_core(self):
        # should be no failing prerequisites
        assert core.check_prereqs(self.config) == []
        # always runs
        assert core.is_enabled(None)

    def test_apply_cluster_rules(self):
        detected_types, cluster_type_hits = hmm_detection.apply_cluster_rules(self.record, self.results_by_id,
                                                                              self.rules)
        for gid in detected_types:
            detected_types[gid] = set(detected_types[gid])
        expected_types = {
            "GENE_1": set(["MetaboliteA", "MetaboliteB", "MetaboliteC", "MetaboliteD"]),
            "GENE_2": set(["MetaboliteC", "MetaboliteD"]),
            "GENE_3": set(["Metabolite0"]),
            "GENE_4": set(["MetaboliteA"]),
            "GENE_5": set(["Metabolite1", "MetaboliteA"])
        }
        assert detected_types == expected_types

        assert cluster_type_hits == {"MetaboliteA": {"GENE_1", "GENE_4", "GENE_5"},
                                     "MetaboliteB": {"GENE_1"},
                                     "MetaboliteC": {"GENE_1", "GENE_2"},
                                     'MetaboliteD': {'GENE_1', 'GENE_2'},
                                     'Metabolite0': {'GENE_3'},
                                     'Metabolite1': {'GENE_5'}}

    def test_find_protoclusters(self):
        cds_features_by_type = {"MetaboliteA": {"GENE_1", "GENE_4", "GENE_5"},
                                "MetaboliteB": {"GENE_1"},
                                "MetaboliteC": {"GENE_1", "GENE_2"},
                                'MetaboliteD': {'GENE_1', 'GENE_2'},
                                'Metabolite0': {'GENE_3'},
                                'Metabolite1': {'GENE_5'}}
        rules = {rule.name: rule for rule in self.rules}
        dummy_extender_domains = {}
        for cluster in hmm_detection.find_protoclusters(self.record, cds_features_by_type, rules,
                                                        self.results_by_id, dummy_extender_domains):
            self.record.add_protocluster(cluster)
        assert len(self.record.get_protoclusters()) == 7
        cluster_products = sorted([cluster.product for cluster in self.record.get_protoclusters()])
        assert cluster_products == sorted([f"Metabolite{i}" for i in "01AABCD"])
        self.record.create_candidate_clusters()
        assert len(self.record.get_candidate_clusters()) == 3
        self.record.create_regions()
        assert len(self.record.get_regions()) == 3
        result_regions = []
        for region in self.record.get_regions():
            result_regions.append(sorted(cds.get_name() for cds in region.cds_children))

        expected_regions = [
            ["GENE_1", "GENE_2"],
            ["GENE_3"],
            ["GENE_4", "GENE_5"]
        ]
        assert result_regions == expected_regions

    def test_create_rules(self):
        aliases = {}
        rules = hmm_detection.create_rules([self.rules_file], self.signature_names,
                                           self.valid_categories, aliases)
        assert len(rules) == open(self.rules_file, encoding="utf-8").read().count("\nRULE")
        t1pks_rules = [rule for rule in rules if rule.name == "T1PKS"]
        assert len(t1pks_rules) == 1
        rule = t1pks_rules[0]
        assert rule.cutoff == 20000
        assert rule.neighbourhood == 20000

    def test_profiles_parsing(self):
        rule_files = [path.get_full_path(__file__, "..", "cluster_rules", f"{name}.txt")
                      for name in ("strict", "relaxed", "loose")]
        profiles_used = set()

        with open(self.filter_file, "r", encoding="utf-8") as handle:
            filter_lines = handle.readlines()
        for line in filter_lines:
            for sig in line.split(','):
                profiles_used.add(sig.strip())

        rules = hmm_detection.create_rules(rule_files, self.signature_names, self.valid_categories)
        for rule in rules:
            profiles_used = profiles_used.union(rule.conditions.profiles)
            for related in rule.related:
                profiles_used.add(related)

        names = self.signature_names

        signatures_not_in_rules = names.difference(profiles_used)
        assert not signatures_not_in_rules, f"No rules use {signatures_not_in_rules}"

        profiles_without_signature = profiles_used.difference(names)
        assert not profiles_without_signature, f"No signature definitions for {profiles_without_signature}"

    def test_filter(self):
        # fake HSPs all in one CDS with overlap > 20 and query_ids from the same equivalence group
        equivalence_groups = [{"AMP-binding", "A-OX"}]

        # not overlapping by > 20
        first = FakeHSPHit("AMP-binding", "A", 50, 90, 0.1, None)
        second = FakeHSPHit("A-OX", "A", 70, 100, 0.5, None)
        new, by_id = hmm_detection.filter_results([first, second], {"A": [first, second]},
                                                  equivalence_groups)
        assert new == [first, second]
        assert by_id == {"A": [first, second]}

        # overlapping, in same group
        first.hit_end = 91
        assert hmm_detection.hsp_overlap_size(first, second) == 21
        new, by_id = hmm_detection.filter_results([first, second], {"A": [first, second]},
                                                  equivalence_groups)
        assert new == [second]
        assert by_id == {"A": [second]}

        # overlapping, not in same group
        second.query_id = "none"
        new, by_id = hmm_detection.filter_results([first, second], {"A": [first, second]},
                                                  equivalence_groups)
        assert new == [first, second]
        assert by_id == {"A": [first, second]}

        # not in the same CDS, but int he same group
        second.hit_id = "B"
        second.query_id = "A-OX"
        new, by_id = hmm_detection.filter_results([first, second], {"A": [first], "B": [second]},
                                                  equivalence_groups)
        assert new == [first, second]
        assert by_id == {"A": [first], "B": [second]}

    def test_filter_multiple(self):
        # all in one CDS no overlap and the same query_ids -> cull all but the best score

        # not overlapping, not same query_id
        first = FakeHSPHit("AMP-binding", "A", 50, 60, 0.1, None)
        second = FakeHSPHit("A-OX", "A", 70, 100, 0.5, None)
        both = [first, second]
        by_id = {"A": [first, second]}
        new, by_id = hmm_detection.filter_result_multiple(list(both), dict(by_id))
        assert new == [first, second]
        assert by_id == {"A": [first, second]}

        # not overlapping, same query_id
        first.query_id = "A-OX"
        new, by_id = hmm_detection.filter_result_multiple(list(both), dict(by_id))
        assert new == [second]
        assert by_id == {"A": [second]}

        # not in same CDS, same query_id
        second.hit_id = "B"
        by_id = {"A": [first], "B": [second]}
        new, by_id = hmm_detection.filter_result_multiple(list(both), dict(by_id))
        assert new == [first, second]
        assert by_id == {"A": [first], "B": [second]}

    def test_equivalence_groups(self):
        group_file = path.get_full_path(os.path.dirname(__file__), "filterhmmdetails.txt")
        sets = []
        with open(group_file, encoding="utf-8") as group_lines:
            sets = [set(line.strip().split(',')) for line in group_lines]

        # ensure they have at least two elements
        assert all(len(s) > 1 for s in sets)

        # ensure that the groups are disjoint
        for i, group in enumerate(sets):
            for other in sets[i + 1:]:
                assert group.isdisjoint(other)

    def test_hsp_overlap_size(self):
        overlap_size = hmm_detection.hsp_overlap_size
        first = FakeHSPHit("A", "A", 50, 60, 0., None)
        second = FakeHSPHit("B", "B", 70, 100, 0., None)
        # no overlap
        assert overlap_size(first, second) == 0
        first.hit_end = 70
        # still no overlap, end isn't inclusive
        assert overlap_size(first, second) == 0
        # a mix of second starting inside first
        for i in range(1, 30):
            first.hit_end += 1
            assert overlap_size(first, second) == i
        # second wholly contained
        first.hit_end = 110
        assert overlap_size(first, second) == 30

        # first inside second
        first.hit_start = 75
        assert overlap_size(first, second) == 25

        # first inside second, but direction reversed
        first.hit_end = 50
        with self.assertRaises(AssertionError):
            overlap_size(first, second)

    def test_hmm_files_and_details_match(self):
        data_dir = path.get_full_path(os.path.dirname(__file__), "data", "")
        details_files = {prof.path for prof in signatures.get_signature_profiles()}
        details_files = {filepath.replace(data_dir, "") for filepath in details_files}
        data_dir_contents = set(glob.glob(data_dir + "*.hmm"))
        data_dir_contents = {filepath.replace(data_dir, "") for filepath in data_dir_contents}
        # ignore bgc_seeds.hmm for the sake of comparison, it's a generated aggregate
        data_dir_contents.discard("bgc_seeds.hmm")
        missing_files = details_files - data_dir_contents
        assert not missing_files
        extra_files = data_dir_contents - details_files
        assert not extra_files
        # finally, just to be sure
        assert data_dir_contents == details_files


class TestRuleExtenders(unittest.TestCase):
    def setUp(self):
        self.rule_name = "MetaboliteA"
        cutoff = 2000
        self.extender_name = "B"
        self.rule_text = "\n".join([
            f"RULE {self.rule_name}",
            "CATEGORY Cat",
            f"CUTOFF {cutoff//1000}",
            "NEIGHBOURHOOD 0",
            f"CONDITIONS A and {self.extender_name}",
            ])
        # two CDSes need to be at the boundaries outside the cutoffs, the rest
        # should be within cutoff of the next CDS
        self.cdses = [
            DummyCDS(locus_tag="X0", start=0, end=500),
            DummyCDS(locus_tag="1", start=3000, end=4000),
            DummyCDS(locus_tag="2", start=5000, end=6000),
            DummyCDS(locus_tag="3", start=7000, end=8000),
            DummyCDS(locus_tag="4", start=9000, end=10000),
            DummyCDS(locus_tag="5", start=11000, end=12000),
            DummyCDS(locus_tag="X1", start=15500, end=16000),
        ]
        self.potential_cores = self.cdses[1:-1]
        assert set("12345") == {cds.get_name() for cds in self.potential_cores}

        self.record = DummyRecord(features=self.cdses)
        assert all(self.record.get_cds_by_name(cds.get_name()) for cds in self.cdses)
        self.results_by_id = {}
        self.add_hit("3", "A")
        self.add_hit("X0", self.extender_name)
        self.add_hit("X1", self.extender_name)

    def detect(self, add_extenders):
        # the X signature/requirement is present just to make a valid CDS condition
        rule_text = self.rule_text
        if add_extenders:
            rule_text += f" EXTENDERS cds({self.extender_name} and not X)"
        rules = rule_parser.Parser(rule_text, {"A", self.extender_name}, {"Cat"}).rules
        assert len(rules) == 1
        if not add_extenders:
            assert not rules[0].extenders
        else:
            assert rules[0].extenders

        pkg = hmm_detection  # to avoid quite a bit of repetition below
        sigs = {name: pkg.HmmSignature(name, "", 0, "") for name in ["A", self.extender_name, "X"]}
        ruleset = create_ruleset(rules, hmm_profiles=sigs)
        with patch.object(pkg, "find_hmmer_hits", return_value=self.results_by_id):
            return pkg.detect_protoclusters_and_signatures(
                self.record, ruleset, "dummy_tool"
            )

    def add_hit(self, cds_name, hit_name):
        hit = FakeHSPHit(hit_name, cds_name, 0, 10, 50, 0)
        if cds_name not in self.results_by_id:
            self.results_by_id[cds_name] = []
        self.results_by_id[cds_name].append(hit)

    def test_no_extension(self):
        # without extenders, 1 and 5 should be too far away from the A in the middle gene
        for name in "1245":
            self.add_hit(name, self.extender_name)
        results = self.detect(add_extenders=False)
        expected = self.potential_cores[1:-1]
        self.check_results(expected, results)

    def test_extenders_left(self):
        # only add extra hits to the left
        for name in [cds.get_name() for cds in self.potential_cores[:2]]:
            self.add_hit(name, self.extender_name)
        expected = self.potential_cores[:3]
        results = self.detect(add_extenders=True)
        self.check_results(expected, results)

    def test_extenders_right(self):
        # only add extra hits to the right
        for name in [cds.get_name() for cds in self.potential_cores[3:]]:
            self.add_hit(name, self.extender_name)

        expected = self.potential_cores[2:]
        results = self.detect(add_extenders=True)
        self.check_results(expected, results)

    def check_results(self, expected, results):
        # the two outsider CDSes should never be in the expected set
        assert self.cdses[0] not in expected
        assert self.cdses[-1] not in expected
        expected = sorted(expected)
        # there should be one protocluster covering all cdses of interest
        assert len(results.protoclusters) == 1
        protocluster = results.protoclusters[0]
        for cds in self.cdses:
            assert cds not in expected or cds.is_contained_by(protocluster)
        assert protocluster.location.start == expected[0].location.start
        assert protocluster.location.end == expected[-1].location.end

        # initially all gene functions should all be OTHER
        for cds in self.cdses:
            assert cds.gene_function == cds.gene_function.OTHER
        results.annotate_cds_features()
        # but once annotated, their gene functions should be CORE even if they were an extender
        for cds in self.cdses:
            if cds in expected:
                assert cds.gene_function == cds.gene_function.CORE
            else:
                assert cds.gene_function == cds.gene_function.OTHER


class TestSignatureFile(unittest.TestCase):
    def test_details(self):
        data_dir = path.get_full_path(os.path.dirname(__file__), 'data')
        check_hmm_signatures(os.path.join(data_dir, 'hmmdetails.txt'), data_dir)


class TestDynamicGather(unittest.TestCase):
    def _go(self, dummy_module):
        with patch.object(pkgutil, "walk_packages", return_value=[Namespace(name="dummy")]):
            with patch.object(importlib, "import_module", return_value=dummy_module):
                return core._get_dynamic_profiles()

    def test_gather(self):
        prof_a = DynamicProfile("A", "desc a", lambda rec: {})
        prof_b = DynamicProfile("b", "desc b", lambda rec: {})
        dynamics = self._go(Namespace(a=prof_a, b=prof_b, c="some text"))
        assert len(dynamics) == 2
        assert sorted(list(dynamics)) == ["A", "b"]
        assert dynamics["A"].name == "A"
        assert dynamics["b"].name == "b"

    def test_duplicates(self):
        prof_a = DynamicProfile("A", "desc a", lambda rec: {})
        prof_b = DynamicProfile("A", "desc b", lambda rec: {})
        with self.assertRaisesRegex(ValueError, "duplicate dynamic profile"):
            self._go(Namespace(a=prof_a, b=prof_b))

    def test_empty(self):
        with self.assertRaisesRegex(ValueError, "subpackage .* has no"):
            self._go(Namespace(a="7"))


class TestRulesetConstruction(unittest.TestCase):
    def get_ruleset(self, names=None, categories=None):
        if names is None:
            names = set()
        if categories is None:
            categories = set()
        args = []
        if names:
            args.extend(["--hmmdetection-limit-to-rule-names", ",".join(names)])
        if categories:
            args.extend(["--hmmdetection-limit-to-rule-categories", ",".join(categories)])
        options = build_config(args, modules=[core])
        assert not core.check_options(options)
        return core.get_ruleset(options)

    def test_no_restrictions(self):
        ruleset = self.get_ruleset()
        assert len(ruleset.rules) > 50  # it's many more, but this is fine

    def test_rules(self):
        subset = {"NRPS", "sactipeptide"}
        ruleset = self.get_ruleset(names=subset)
        assert len(ruleset.rules) == len(subset)
        for rule in ruleset.rules:
            assert rule.name in subset

    def test_categories(self):
        subset = {"terpene", "PKS"}
        ruleset = self.get_ruleset(categories=subset)
        assert len(ruleset.rules) > 5  # again, the exact number isn't important
        for rule in ruleset.rules:
            assert rule.category in subset

    def test_intersection(self):
        categories = {"NRPS"}
        rules = {"NRPS", "sactipeptide"}
        ruleset = self.get_ruleset(names=rules, categories=categories)
        assert len(ruleset.rules) == 1
        for rule in ruleset.rules:
            assert rule.category in categories
            assert rule.name in rules


class TestOptions(unittest.TestCase):
    def test_invalid_rule_names(self):
        config = build_config([
            "--hmmdetection-limit-to-rule-names", ",".join({"NRPS", "xyz_rule"}),
        ], modules=[core])
        issues = core.check_options(config)
        assert len(issues) == 1
        assert issues[0] == "Unknown rules in requested rule subset: {'xyz_rule'}"

    def test_invalid_rule_categories(self):
        config = build_config([
            "--hmmdetection-limit-to-rule-categories", ",".join({"NRPS", "xyz_cat"}),
        ], modules=[core])
        issues = core.check_options(config)
        assert len(issues) == 1
        assert issues[0] == "Unknown rules in requested rule category subset: {'xyz_cat'}"

    def test_combined(self):
        config = build_config([
            "--hmmdetection-limit-to-rule-names", ",".join({"NRPS", "xyz_rule"}),
            "--hmmdetection-limit-to-rule-categories", ",".join({"NRPS", "xyz_cat"}),
        ], modules=[core])
        issues = core.check_options(config)
        assert len(issues) == 2
        assert set(issues) == {
            "Unknown rules in requested rule subset: {'xyz_rule'}",
            "Unknown rules in requested rule category subset: {'xyz_cat'}",
        }


class TestMultipliers(unittest.TestCase):
    def setUp(self):
        self.record = DummyRecord()
        self.record.add_cds_feature(DummyCDS())
        rule_text = "RULE MetaboliteA CATEGORY Cat CUTOFF 10 NEIGHBOURHOOD 5 CONDITIONS modelA"
        self.signatures = {name: hmm_detection.HmmSignature(name, "", 0, "") for name in ["modelA"]}
        self.rules = rule_parser.Parser(rule_text, set(self.signatures), {"Cat"}).rules

    def tearDown(self):
        destroy_config()

    def run_through(self, options):
        ruleset = create_ruleset(self.rules, hmm_profiles=self.signatures)
        with patch.object(core, "detect_protoclusters_and_signatures",
                          side_effect=RuntimeError("stop here")) as patched:
            with patch.object(core, "get_rule_categories", return_value=[]):
                with self.assertRaisesRegex(RuntimeError, "stop here"):
                    core.run_on_record(self.record, None, options)
            assert patched.called_once
            # find the ruleset arg used, and if it doesn't exist, failing here is fine
            ruleset = [arg for arg in patched.call_args[0] if isinstance(arg, hmm_detection.Ruleset)][0]
        return ruleset

    def test_bacteria_ignores_fungal_multi(self):
        options = build_config([
            "--taxon", "bacteria",
            "--hmmdetection-fungal-cutoff-multiplier", "7.0",
            "--hmmdetection-fungal-neighbourhood-multiplier", "3",
        ], modules=[core])
        ruleset = self.run_through(options)
        assert ruleset.multipliers == core.Multipliers(1.0, 1.0)

    def test_fungal_respects_multis(self):
        options = build_config([
            "--taxon", "fungi",
            "--hmmdetection-fungal-cutoff-multiplier", "0.5",
            "--hmmdetection-fungal-neighbourhood-multiplier", "3",
        ], modules=[core])
        ruleset = self.run_through(options)
        assert ruleset.multipliers == core.Multipliers(0.5, 3.0)

    def test_check_options(self):
        options = build_config([
            "--taxon", "fungi",
            "--hmmdetection-fungal-cutoff-multiplier", "0",
            "--hmmdetection-fungal-neighbourhood-multiplier", "-5",
        ], modules=[core])
        assert core.check_options(options) == [
            "Invalid fungal cutoff multiplier: 0.0",
            "Invalid fungal neighbourhood multiplier: -5.0",
        ]

    def test_reuse_changes(self):
        options = build_config([
            "--taxon", "fungi",
            "--hmmdetection-fungal-cutoff-multiplier", "1",
            "--hmmdetection-fungal-neighbourhood-multiplier", "1.5",
        ], modules=[core])
        results = core.run_on_record(self.record, None, options)
        as_json = json.loads(json.dumps(results.to_json()))

        regenerated = core.regenerate_previous_results(as_json, self.record, options)
        assert regenerated.rule_results.multipliers == results.rule_results.multipliers

        # ensure a changed cutoff multiplier breaks results reuse
        options = build_config([
            "--taxon", "fungi",
            "--hmmdetection-fungal-cutoff-multiplier", "2.0",
        ], modules=[core])
        with self.assertRaisesRegex(RuntimeError, "cutoff multiplier .* incompatible"):
            core.regenerate_previous_results(as_json, self.record, options)

        options = build_config([
            "--taxon", "fungi",
            "--hmmdetection-fungal-neighbourhood-multiplier", "0.5",
        ], modules=[core])
        with self.assertRaisesRegex(RuntimeError, "neighbourhood multiplier .* incompatible"):
            core.regenerate_previous_results(as_json, self.record, options)
