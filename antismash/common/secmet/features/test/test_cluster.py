# License: GNU Affero General Public License v3 or later
# A copy of GNU AGPL v3 should have been included in this software package in LICENSE.txt.

# for test files, silence irrelevant and noisy pylint warnings
# pylint: disable=no-self-use,protected-access,missing-docstring

import unittest

from antismash.common.secmet import FeatureLocation
from antismash.common.secmet.features import Cluster


def create_cluster():
    cluster = Cluster(FeatureLocation(8, 71, strand=1),
                      FeatureLocation(3, 76, strand=1), tool="test",
                      cutoff=17, neighbourhood_range=5, product='a')
    cluster.detection_rule = "some rule text"
    return cluster


class TestCluster(unittest.TestCase):
    def setUp(self):
        self.cluster = create_cluster()

    def test_orphaned_numbering(self):
        with self.assertRaisesRegex(ValueError, "Cluster not in a record"):
            print(self.cluster.get_cluster_number())

    def test_biopython_conversion(self):
        bio = self.cluster.to_biopython()
        assert len(bio) == 2
        assert bio[0].type == "cluster" and bio[1].type == "cluster_core"
        new = Cluster.from_biopython(bio[0])
        assert new is not self.cluster
        assert new.cutoff == self.cluster.cutoff == 17
        assert new.neighbourhood_range == self.cluster.neighbourhood_range == 5
        assert new.product == self.cluster.product == 'a'
        assert new.location.start == self.cluster.location.start == 3
        assert new.core_location.start == self.cluster.core_location.start == 8
        assert new.detection_rule == self.cluster.detection_rule == "some rule text"
        assert new.tool == self.cluster.tool == "test"
        assert new.location.start == self.cluster.location.start == 3
        assert new.location.end == self.cluster.location.end == 76
        assert new.core_location.start == self.cluster.core_location.start == 8
        assert new.core_location.end == self.cluster.core_location.end == 71
