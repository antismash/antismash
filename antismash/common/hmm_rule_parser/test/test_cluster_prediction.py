# License: GNU Affero General Public License v3 or later
# A copy of GNU AGPL v3 should have been included in this software package in LICENSE.txt.

# for test files, silence irrelevant and noisy pylint warnings
# pylint: disable=no-self-use,protected-access,missing-docstring,too-many-public-methods

import unittest

from antismash.common.hmm_rule_parser import cluster_prediction, rule_parser
from antismash.common.secmet.features import Cluster, FeatureLocation


class DummyConditions(rule_parser.Conditions):
    """ so a DetectionRule can be created without failing its internal checks """
    def __init__(self):
        super().__init__(negated=False)

    def contains_positive_condition(self):
        return True

# NOTE: the rest of the cluster_prediction tests are still in hmm_detection tests


class TestRedundancy(unittest.TestCase):
    def setUp(self):
        superior = rule_parser.DetectionRule("superior", 10, 10, DummyConditions())
        inferior = rule_parser.DetectionRule("inferior", 10, 10, DummyConditions(), superiors=["superior"])
        irrelevant = rule_parser.DetectionRule("irrelevant", 10, 10, DummyConditions())
        self.rules_by_name = {rule.name: rule for rule in [superior, inferior, irrelevant]}

    def remove(self, clusters):
        return cluster_prediction.remove_redundant_clusters(clusters, self.rules_by_name)

    def create_cluster(self, rule_name, start, end):
        rule = self.rules_by_name[rule_name]
        core = FeatureLocation(start, end)
        surrounds = FeatureLocation(max(0, start - rule.extent), end + rule.extent)
        return Cluster(core, surrounds, tool="testing", cutoff=rule.cutoff,
                       neighbourhood_range=rule.extent, product=rule_name,
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
        print("testing clusters:", clusters)
        assert self.remove(clusters) == clusters

    def test_extents_dont_matter(self):
        extent = self.rules_by_name["superior"].extent
        for new_extent in [extent - 10, extent + 10]:
            self.rules_by_name["inferior"].extent = new_extent
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
