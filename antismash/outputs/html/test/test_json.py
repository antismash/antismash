# License: GNU Affero General Public License v3 or later
# A copy of GNU AGPL v3 should have been included in this software package in LICENSE.txt.

# for test files, silence irrelevant and noisy pylint warnings
# pylint: disable=no-self-use,protected-access,missing-docstring

import unittest

from antismash.common.test.helpers import DummyCluster
from antismash.outputs.html import js as jsonify

class TestClusterPacking(unittest.TestCase):
    def check(self, clusters, expected, padding=None):
        if padding is not None:
            assert jsonify._find_non_overlapping_cluster_groups(clusters, padding) == expected
        else:
            assert jsonify._find_non_overlapping_cluster_groups(clusters) == expected

    def test_no_clusters(self):
        self.check([], {})

    def test_single(self):
        cluster = DummyCluster()
        self.check([cluster], {cluster: 0})

    def test_overlap(self):
        clusters = [DummyCluster(start=100, end=500), DummyCluster(start=300, end=700)]
        self.check(clusters, {clusters[0]: 0, clusters[1]: 1})

    def test_non_overlap(self):
        clusters = [DummyCluster(start=100, end=500), DummyCluster(start=800, end=900)]
        self.check(clusters, {clusters[0]: 0, clusters[1]: 0})

    def test_padding(self):
        # default of 100
        clusters = [DummyCluster(start=100, end=500), DummyCluster(start=550, end=900)]
        self.check(clusters, {clusters[0]: 0, clusters[1]: 1})
        self.check(clusters, {clusters[0]: 0, clusters[1]: 1}, padding=100)

        self.check(clusters, {clusters[0]: 0, clusters[1]: 0}, padding=49)
        self.check(clusters, {clusters[0]: 0, clusters[1]: 1}, padding=50)
        self.check(clusters, {clusters[0]: 0, clusters[1]: 1}, padding=51)

    def test_mixed(self):
        clusters = [
            DummyCluster(start=100, end=190),
            DummyCluster(start=150, end=800),
            DummyCluster(start=200, end=300),
            DummyCluster(start=250, end=500),
            DummyCluster(start=400, end=800),
            DummyCluster(start=750, end=900),
        ]
        expected = {
            clusters[0]: 0,
            clusters[1]: 1,
            clusters[2]: 0,
            clusters[3]: 2,
            clusters[4]: 0,
            clusters[5]: 2,
        }
        self.check(clusters, expected, padding=1)

    def test_bad_padding(self):
        with self.assertRaisesRegex(ValueError, "cannot be negative"):
            self.check([], {}, padding=-1)
