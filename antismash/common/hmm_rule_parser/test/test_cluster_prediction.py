# License: GNU Affero General Public License v3 or later
# A copy of GNU AGPL v3 should have been included in this software package in LICENSE.txt.

# for test files, silence irrelevant and noisy pylint warnings
# pylint: disable=use-implicit-booleaness-not-comparison,protected-access,missing-docstring,too-many-public-methods

import unittest
from unittest.mock import patch

from antismash.common.hmm_rule_parser import cluster_prediction, rule_parser, structures
from antismash.common.secmet.features import Protocluster, FeatureLocation
from antismash.common.test.helpers import DummyRecord, DummyCDS


class DummyConditions(rule_parser.Conditions):
    """ so a DetectionRule can be created without failing its internal checks """
    def __init__(self):
        super().__init__(negated=False)

    def contains_positive_condition(self):
        return True

# NOTE: the rest of the cluster_prediction tests are still in hmm_detection tests


class TestRedundancy(unittest.TestCase):
    def setUp(self):
        superior = rule_parser.DetectionRule("superior", "category", 10, 10, DummyConditions())
        inferior = rule_parser.DetectionRule("inferior", "category", 10, 10, DummyConditions(), superiors=["superior"])
        irrelevant = rule_parser.DetectionRule("irrelevant", "category", 10, 10, DummyConditions())
        self.rules_by_name = {rule.name: rule for rule in [superior, inferior, irrelevant]}

    def remove(self, clusters):
        return cluster_prediction.remove_redundant_protoclusters(clusters, self.rules_by_name)

    def create_cluster(self, rule_name, start, end):
        rule = self.rules_by_name[rule_name]
        core = FeatureLocation(start, end)
        surrounds = FeatureLocation(max(0, start - rule.neighbourhood), end + rule.neighbourhood)
        return Protocluster(core, surrounds, tool="testing", cutoff=rule.cutoff,
                            neighbourhood_range=rule.neighbourhood, product=rule_name,
                            detection_rule="rule text")

    def test_alone(self):
        clusters = [self.create_cluster("inferior", 101, 110)]
        assert clusters == self.remove(clusters)

    def test_non_overlap(self):
        clusters = [self.create_cluster("inferior", 101, 110),
                    self.create_cluster("superior", 200, 210)]
        assert clusters == self.remove(clusters)

    def test_not_relevant_equal(self):
        clusters = [self.create_cluster("inferior", 101, 110),
                    self.create_cluster("irrelevant", 101, 110)]
        assert clusters == self.remove(clusters)

    def test_not_relevant_contained(self):
        clusters = [self.create_cluster("inferior", 105, 108),
                    self.create_cluster("irrelevant", 101, 110)]
        assert clusters == self.remove(clusters)

    def test_not_relevant_larger(self):
        clusters = [self.create_cluster("inferior", 101, 110),
                    self.create_cluster("irrelevant", 105, 108)]
        assert clusters == self.remove(clusters)

    def test_contained(self):
        clusters = [self.create_cluster("inferior", 110, 210),
                    self.create_cluster("superior", 101, 310)]
        assert self.remove(clusters) == [clusters[1]]

    def test_equal(self):
        clusters = [self.create_cluster("inferior", 101, 110),
                    self.create_cluster("superior", 101, 110)]
        assert self.remove(clusters) == [clusters[1]]

    def test_larger(self):
        clusters = [self.create_cluster("inferior", 101, 110),
                    self.create_cluster("superior", 102, 109)]
        assert self.remove(clusters) == [clusters[1]]

    def test_neighbourhoods_dont_matter(self):
        neighbourhood = self.rules_by_name["superior"].neighbourhood
        for new_neighbourhood in [neighbourhood - 10, neighbourhood + 10]:
            self.rules_by_name["inferior"].neighbourhood = new_neighbourhood
            self.test_larger()
            self.test_equal()
            self.test_contained()

    def test_cutoffs_dont_matter(self):
        cutoff = self.rules_by_name["superior"].cutoff
        for new_cutoff in [cutoff - 10, cutoff + 10]:
            self.rules_by_name["inferior"].cutoff = new_cutoff
            self.test_larger()
            self.test_equal()
            self.test_contained()


class TestDynamic(unittest.TestCase):
    def test_find_dynamic(self):
        expected_a = {"cds_name": [structures.DynamicHit("prof_a", "cds_name")]}
        expected_b = {"cds_name": [structures.DynamicHit("prof_b", "cds_name")]}
        profile_a = structures.DynamicProfile("prof_a", "desc a", lambda record: expected_a)
        profile_b = structures.DynamicProfile("prof_b", "desc b", lambda record: expected_b)
        results = cluster_prediction.find_dynamic_hits(DummyRecord(), [profile_a, profile_b])
        assert results["cds_name"] == expected_a["cds_name"] + expected_b["cds_name"]

    @patch.object(cluster_prediction, "find_hmmer_hits", return_value={})
    @patch.object(cluster_prediction, "get_signature_profiles", return_value={})
    def test_full(self, _patched_find, _patched_sigs):
        cdses = [DummyCDS(locus_tag="A", start=0, end=6), DummyCDS(locus_tag="B", start=6, end=12)]
        record = DummyRecord(features=cdses)

        # create a dummy dynamic profile
        def find_a(rec):
            hits = {}
            for cds in rec.get_cds_features():
                if cds.get_name() == "A":
                    hits[cds.get_name()] = [structures.DynamicHit(cds.get_name(), "a_finder")]
            return hits
        profile = structures.DynamicProfile("a_finder", "desc", find_a)
        # make sure the 'profile' functions as expected
        assert cdses[0].get_name() in profile.find_hits(record)

        # build a dummy rule that will search for this hit
        condition = rule_parser.SingleCondition(False, "a_finder")
        rule = rule_parser.DetectionRule("test-name", "Other", 5000, 5000, condition)
        with patch.object(cluster_prediction, "create_rules", return_value=[rule]):
            results = cluster_prediction.detect_protoclusters_and_signatures(
                record, None, None, [None], set("Other"), None, "test_tool",
                dynamic_profiles={profile.name: profile}
            )
        assert results
        assert results.cds_by_cluster
        assert results.protoclusters
        proto = results.protoclusters[0]
        assert proto.product == "test-name"

        results.annotate_cds_features()
        assert cdses[0].sec_met.domains[0].name == "a_finder"

    def test_overlap_names(self):
        record = DummyRecord(features=[DummyCDS()])
        profile = structures.DynamicProfile("dummy", "desc", lambda record: {})
        with patch.object(cluster_prediction, "get_signature_profiles", return_value=[profile]):
            with self.assertRaisesRegex(ValueError, "profiles overlap"):
                cluster_prediction.detect_protoclusters_and_signatures(
                    record, None, None, [None], set("Other"), None, "test_tool",
                    dynamic_profiles={profile.name: profile}
                )


class TestMultipliers(unittest.TestCase):
    def test_create_rules(self):
        text = "RULE A CATEGORY Cat CUTOFF 10 NEIGHBOURHOOD 5 CONDITIONS A"
        # with default multipliers
        with patch("builtins.open", unittest.mock.mock_open(read_data=text)):
            rule = cluster_prediction.create_rules("dummy.file", {"A"}, {"Cat"}, {})[0]
        assert rule.cutoff == 10_000
        assert rule.neighbourhood == 5_000

        # with custom multipliers
        multipliers = cluster_prediction.Multipliers(
            cutoff=1.5,
            neighbourhood=0.5,
        )
        with patch("builtins.open", unittest.mock.mock_open(read_data=text)):
            multiplied = cluster_prediction.create_rules("dummy.file", {"A"}, {"Cat"}, {},
                                                         multipliers=multipliers,
                                                         )[0]
        # make sure the multipliers were used
        assert multiplied.cutoff == rule.cutoff * multipliers.cutoff
        assert multiplied.neighbourhood == rule.neighbourhood * multipliers.neighbourhood

    @patch.object(cluster_prediction, "get_signature_profiles", return_value={})
    def test_multplier_propagation(self, _patched_sig):
        record = DummyRecord()
        record.add_cds_feature(DummyCDS())
        args = [record, "dummy.sigs", "dummy.seeds", ["dummy.rules"], {"cat"},
                "dummy.filter", "tool"]
        multipliers = cluster_prediction.Multipliers(
            cutoff=0.1,
            neighbourhood=2.0,
        )
        kwargs = {
            "multipliers": multipliers,
        }
        with patch.object(cluster_prediction, "create_rules",
                          side_effect=RuntimeError("stop here")) as patched_create:
            with self.assertRaisesRegex(RuntimeError, "stop here"):
                cluster_prediction.detect_protoclusters_and_signatures(*args, **kwargs)
            # make sure the multipliers made it all the way to rule creation
            _, actual_kwargs = patched_create.call_args
            assert actual_kwargs["multipliers"] == multipliers
