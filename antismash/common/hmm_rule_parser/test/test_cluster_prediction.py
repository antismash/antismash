# License: GNU Affero General Public License v3 or later
# A copy of GNU AGPL v3 should have been included in this software package in LICENSE.txt.

# for test files, silence irrelevant and noisy pylint warnings
# pylint: disable=no-self-use,protected-access,missing-docstring,too-many-public-methods

import unittest

from antismash.common.hmm_rule_parser import cluster_prediction, rule_parser
from antismash.common.secmet.features import ClusterBorder, FeatureLocation


class DummyConditions(rule_parser.Conditions):
    """ so a DetectionRule can be created without failing its internal checks """
    def __init__(self):
        super().__init__(negated=False)

    def contains_positive_condition(self):
        return True

# NOTE: the rest of the cluster_prediction tests are still in hmm_detection tests


class TestRedundancy(unittest.TestCase):
    def setUp(self):
        superior = rule_parser.DetectionRule("superior", 10000, 10000, DummyConditions())
        inferior = rule_parser.DetectionRule("inferior", 10000, 10000, DummyConditions(), superiors=["superior"])
        irrelevant = rule_parser.DetectionRule("irrelevant", 10000, 10000, DummyConditions())
        self.rules_by_name = {rule.name: rule for rule in [superior, inferior, irrelevant]}

    def remove(self, borders):
        return cluster_prediction.remove_redundant_borders(borders, self.rules_by_name)

    def create_border(self, rule_name, start, end):
        rule = self.rules_by_name[rule_name]
        return ClusterBorder(FeatureLocation(start, end), tool="testing",
                             cutoff=rule.cutoff, extent=rule.extent, product=rule_name)

    def test_alone(self):
        borders = [self.create_border("inferior", 1, 10)]
        assert borders == self.remove(borders)

    def test_non_overlap(self):
        borders = [self.create_border("inferior", 1, 10),
                   self.create_border("superior", 100, 110)]
        assert borders == self.remove(borders)

    def test_not_relevant_equal(self):
        borders = [self.create_border("inferior", 1, 10),
                   self.create_border("irrelevant", 1, 10)]
        assert borders == self.remove(borders)

    def test_not_relevant_contained(self):
        borders = [self.create_border("inferior", 5, 8),
                   self.create_border("irrelevant", 1, 10)]
        assert borders == self.remove(borders)

    def test_not_relevant_larger(self):
        borders = [self.create_border("inferior", 1, 10),
                   self.create_border("irrelevant", 5, 8)]
        assert borders == self.remove(borders)

    def test_contained(self):
        borders = [self.create_border("inferior", 10, 110),
                   self.create_border("superior", 1, 210)]
        assert self.remove(borders) == [borders[1]]

    def test_equal(self):
        borders = [self.create_border("inferior", 1, 10),
                   self.create_border("superior", 1, 10)]
        assert self.remove(borders) == [borders[1]]

    def test_larger(self):
        borders = [self.create_border("inferior", 1, 10),
                   self.create_border("superior", 2, 9)]
        assert self.remove(borders) == borders

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
