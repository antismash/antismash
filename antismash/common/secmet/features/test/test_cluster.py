# License: GNU Affero General Public License v3 or later
# A copy of GNU AGPL v3 should have been included in this software package in LICENSE.txt.

# for test files, silence irrelevant and noisy pylint warnings
# pylint: disable=no-self-use,protected-access,missing-docstring

import unittest

from antismash.common.secmet import FeatureLocation
from antismash.common.secmet.features import Cluster, CDSFeature
from antismash.common.secmet.qualifiers import SecMetQualifier


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


class TestDefinitionCDS(unittest.TestCase):
    def setUp(self):
        self.cluster = create_cluster()
        self.cluster.core_location = FeatureLocation(30, 50)
        self.inside_cds = CDSFeature(FeatureLocation(40, 45, strand=1), locus_tag="test1")
        self.neighbour_cds = CDSFeature(FeatureLocation(20, 25, strand=1), locus_tag="test2")
        self.outside_cds = CDSFeature(FeatureLocation(120, 125, strand=1), locus_tag="test3")
        assert not self.cluster.cds_children
        assert not self.cluster.definition_cdses

        self.matching_qual = SecMetQualifier(set(["a", "b"]), [SecMetQualifier.Domain("a-b", 1e-5, 1., 2, "test")])
        self.mismatch_qual = SecMetQualifier(set(["b", "c"]), [SecMetQualifier.Domain("b-c", 1e-5, 1., 2, "test")])

    def test_no_secmet_qual_inside(self):
        self.cluster.add_cds(self.inside_cds)
        assert self.cluster.cds_children
        assert not self.cluster.definition_cdses

    def test_matching_secmet_inside(self):
        self.inside_cds.sec_met = self.matching_qual
        self.cluster.add_cds(self.inside_cds)
        assert self.cluster.cds_children
        assert self.cluster.definition_cdses == {self.inside_cds}

    def test_matching_secmet_neighbour(self):
        self.neighbour_cds.sec_met = self.matching_qual
        self.cluster.add_cds(self.neighbour_cds)
        assert self.cluster.cds_children
        assert not self.cluster.definition_cdses

    def test_matching_secmet_outside(self):
        self.outside_cds.sec_met = self.matching_qual
        with self.assertRaisesRegex(ValueError, "not contained by"):
            self.cluster.add_cds(self.outside_cds)

    def test_mismatching_secmet_inside(self):
        self.inside_cds.sec_met = self.mismatch_qual
        self.cluster.add_cds(self.inside_cds)
        assert self.cluster.cds_children
        assert not self.cluster.definition_cdses

    def test_mismatching_secmet_neighbour(self):
        self.neighbour_cds.sec_met = self.mismatch_qual
        self.cluster.add_cds(self.neighbour_cds)
        assert self.cluster.cds_children
        assert not self.cluster.definition_cdses

    def test_mismatching_secmet_outside(self):
        self.outside_cds.sec_met = self.mismatch_qual
        with self.assertRaisesRegex(ValueError, "not contained by"):
            self.cluster.add_cds(self.outside_cds)
