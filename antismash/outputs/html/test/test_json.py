# License: GNU Affero General Public License v3 or later
# A copy of GNU AGPL v3 should have been included in this software package in LICENSE.txt.

# for test files, silence irrelevant and noisy pylint warnings
# pylint: disable=use-implicit-booleaness-not-comparison,protected-access,missing-docstring

import unittest

from antismash.common.test.helpers import (
    DummyCandidateCluster,
    DummyProtocluster,
    DummyRecord,
    DummyRegion,
)
from antismash.outputs.html import js as jsonify


class TestClusterPacking(unittest.TestCase):
    def check(self, clusters, expected, padding=None, record_length=None):
        if padding is not None:
            groups = jsonify._find_non_overlapping_cluster_groups(clusters, padding, record_length=record_length)
        else:
            groups = jsonify._find_non_overlapping_cluster_groups(clusters, record_length=record_length)
        assert groups == expected

    def test_no_clusters(self):
        self.check([], {})

    def test_single(self):
        cluster = DummyProtocluster()
        self.check([cluster], {cluster: 0})

    def test_overlap(self):
        clusters = [DummyProtocluster(start=100, end=500), DummyProtocluster(start=300, end=700)]
        self.check(clusters, {clusters[0]: 0, clusters[1]: 1})

    def test_non_overlap(self):
        clusters = [DummyProtocluster(start=100, end=500), DummyProtocluster(start=800, end=900)]
        self.check(clusters, {clusters[0]: 0, clusters[1]: 0})

    def test_padding(self):
        # default of 100
        clusters = [DummyProtocluster(start=100, end=500), DummyProtocluster(start=550, end=900)]
        self.check(clusters, {clusters[0]: 0, clusters[1]: 1})
        self.check(clusters, {clusters[0]: 0, clusters[1]: 1}, padding=100)

        self.check(clusters, {clusters[0]: 0, clusters[1]: 0}, padding=49)
        self.check(clusters, {clusters[0]: 0, clusters[1]: 1}, padding=50)
        self.check(clusters, {clusters[0]: 0, clusters[1]: 1}, padding=51)

    def test_mixed(self):
        clusters = [
            DummyProtocluster(start=100, end=190),
            DummyProtocluster(start=150, end=800),
            DummyProtocluster(start=200, end=300),
            DummyProtocluster(start=250, end=500),
            DummyProtocluster(start=400, end=800),
            DummyProtocluster(start=750, end=900),
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

    def test_with_cross_origin(self):
        clusters = [
            DummyProtocluster(start=59, end=6, core_start=79, core_end=80, record_length=100),
            DummyProtocluster(start=70, end=39, core_start=90, core_end=19, record_length=100),
            DummyProtocluster(start=70, end=28, core_start=90, core_end=8, record_length=100),
            DummyProtocluster(start=16, end=38, core_start=26, core_end=28, record_length=100),
        ]
        expected = {
            clusters[0]: 0,
            clusters[1]: 1,
            clusters[2]: 2,
            clusters[3]: 0,
        }
        self.check(clusters, expected, record_length=100)
        self.check(sorted(clusters), expected, record_length=100)


class TestConversion(unittest.TestCase):
    def test_neighbourhood_bounded(self):
        length = 100
        record = DummyRecord(seq="A" * length)
        padding = length // 4
        protocluster = DummyProtocluster(start=0, core_start=padding, core_end=length - padding,
                                         end=length, neighbourhood_range=padding*2)
        record.add_protocluster(protocluster)
        # the protocluster must fit, so that the conversion is tested with decent values
        assert protocluster.location.end == length
        # and a naive core end + neighbourhood size must be outside the record, to test the conversion
        assert (protocluster.core_location.end + protocluster.neighbourhood_range) > length
        assert len(protocluster.location) == length
        candidate = DummyCandidateCluster([protocluster])
        record.add_candidate_cluster(candidate)
        region = DummyRegion(candidate_clusters=[candidate])
        record.add_region(region)

        assert not record.get_subregions()

        # first, try non-circular records, checking that they're bounded correctly
        result = jsonify.get_clusters_from_region(region, record_length=100, circular=False)
        types = [area["kind"] for area in result]
        assert types == ["candidatecluster", "subregion", "protocluster"]
        result = result[-1]
        assert result["neighbouring_start"] == 0
        assert result["neighbouring_end"] == length

        # then, try circular records, checking that they can extend past the origin
        result = jsonify.get_clusters_from_region(region, record_length=100, circular=True)
        types = [area["kind"] for area in result]
        assert types == ["candidatecluster", "subregion", "protocluster"]
        result = result[-1]
        assert result["neighbouring_end"] == protocluster.core_location.end + protocluster.neighbourhood_range
        assert result["neighbouring_start"] == -padding
