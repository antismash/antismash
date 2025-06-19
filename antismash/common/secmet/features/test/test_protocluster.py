# License: GNU Affero General Public License v3 or later
# A copy of GNU AGPL v3 should have been included in this software package in LICENSE.txt.

# for test files, silence irrelevant and noisy pylint warnings
# pylint: disable=use-implicit-booleaness-not-comparison,protected-access,missing-docstring

import unittest

from antismash.common.secmet import FeatureLocation
from antismash.common.secmet.features.protocluster import Protocluster, SideloadedProtocluster
from antismash.common.secmet.qualifiers import GeneFunction
from antismash.common.secmet.test.helpers import DummyCDS, DummyRecord


def create_cluster():
    cluster = Protocluster(FeatureLocation(8, 71, strand=1),
                           FeatureLocation(3, 76, strand=1), tool="test",
                           cutoff=17, neighbourhood_range=5, product='a',
                           detection_rule="some rule text",
                           product_category="some category")
    return cluster


class TestProtocluster(unittest.TestCase):
    def setUp(self):
        self.cluster = create_cluster()

    def test_orphaned_numbering(self):
        with self.assertRaisesRegex(ValueError, "Protocluster not in a record"):
            print(self.cluster.get_protocluster_number())

    def test_biopython_conversion(self):
        bio = self.cluster.to_biopython()
        assert len(bio) == 2
        assert bio[0].type == "protocluster" and bio[1].type == "proto_core"
        new = Protocluster.from_biopython(bio[0])
        assert new is not self.cluster
        assert new.cutoff == self.cluster.cutoff == 17
        assert new.neighbourhood_range == self.cluster.neighbourhood_range == 5
        assert new.product == self.cluster.product == 'a'
        assert new.location.start == self.cluster.location.start == 3
        assert new.core_location.start == self.cluster.core_location.start == 8
        assert new.detection_rule == self.cluster.detection_rule == "some rule text"
        assert new.product_category == self.cluster.product_category == "some category"
        assert new.tool == self.cluster.tool == "test"
        assert new.location.start == self.cluster.location.start == 3
        assert new.location.end == self.cluster.location.end == 76
        assert new.core_location.start == self.cluster.core_location.start == 8
        assert new.core_location.end == self.cluster.core_location.end == 71

    def test_contig_edge_start(self):
        cluster = Protocluster(FeatureLocation(20, 30, strand=1),
                               FeatureLocation(10, 40, strand=1), tool="test",
                               cutoff=3, neighbourhood_range=5, product='a',
                               detection_rule="some rule text",
                               product_category="some category")
        # check default
        record = DummyRecord(seq="A" * 100)
        cluster.parent_record = record
        assert cluster.contig_edge is False

        # make the cutoff large enough to reach the edge
        cluster.cutoff = 25
        assert cluster.core_location.start - cluster.cutoff < 0 < cluster.location.start
        assert cluster.contig_edge is True

    def test_contig_edge_end(self):
        cluster = Protocluster(FeatureLocation(50, 60, strand=1),
                               FeatureLocation(40, 70, strand=1), tool="test",
                               cutoff=3, neighbourhood_range=5, product='a',
                               detection_rule="some rule text",
                               product_category="some category")
        # check default
        record = DummyRecord(seq="A" * 80)
        cluster.parent_record = record
        assert cluster.contig_edge is False

        # make the cutoff large enough to reach the edge
        cluster.cutoff = 25
        assert cluster.location.start < len(record.seq) < cluster.core_location.end + cluster.cutoff
        assert cluster.contig_edge is True

    def test_sorting(self):
        location = FeatureLocation(1, 4, 1)

        def create(loc, product):
            return Protocluster(loc, loc, tool="test", cutoff=3, neighbourhood_range=5,
                                product=product, detection_rule="text", product_category="cat")

        first = create(location, "A")
        second = create(location, "B")
        different = create(FeatureLocation(3, 5, 1), "A")
        assert first < second
        assert second < different
        assert second > first


class TestSideloaded(unittest.TestCase):
    def test_conversion(self):
        core = FeatureLocation(8, 71, strand=1)
        surrounds = FeatureLocation(3, 76, strand=1)
        extras = {"a": ["5", "c"], "b": ["something"]}
        source = SideloadedProtocluster(core, surrounds, "tool name", "some-product", extra_qualifiers=extras)

        assert source.neighbourhood_range == 5

        bio_features = source.to_biopython()
        assert len(bio_features) == 2
        for key, val in extras.items():
            assert bio_features[0].qualifiers[key] == val

        for regenerator in [SideloadedProtocluster, Protocluster]:
            dest = regenerator.from_biopython(bio_features[0])
            assert isinstance(dest, SideloadedProtocluster)
            assert dest.extra_qualifiers == source.extra_qualifiers == extras
            assert dest.tool == source.tool
            assert dest.product == source.product
            assert dest.location == source.location
            assert dest.core_location == source.core_location
            assert dest.neighbourhood_range == source.neighbourhood_range

            for key, val in extras.items():
                assert not dest.get_qualifier(key)


class TestDefinitionCDS(unittest.TestCase):
    def setUp(self):
        self.cluster = create_cluster()
        self.cluster._core_location = FeatureLocation(30, 50)
        self.inside_cds = DummyCDS(40, 45)
        self.neighbour_cds = DummyCDS(20, 25)
        self.outside_cds = DummyCDS(120, 125)
        assert not self.cluster.cds_children
        assert not self.cluster.definition_cdses

    def add_core_function(self, cds, cluster_product=True):
        product = self.cluster.product if cluster_product else f"not{self.cluster.product}"
        cds.gene_functions.add(GeneFunction.CORE, "test", "dummy", product)

    def test_no_secmet_qual_inside(self):
        self.cluster.add_cds(self.inside_cds)
        assert self.cluster.cds_children
        assert not self.cluster.definition_cdses

    def test_matching_secmet_inside(self):
        self.add_core_function(self.inside_cds)
        self.cluster.add_cds(self.inside_cds)
        assert self.cluster.cds_children
        assert self.cluster.definition_cdses == {self.inside_cds}

    def test_matching_secmet_neighbour(self):
        self.add_core_function(self.neighbour_cds)
        self.cluster.add_cds(self.neighbour_cds)
        assert self.cluster.cds_children
        assert not self.cluster.definition_cdses

    def test_matching_secmet_outside(self):
        self.add_core_function(self.outside_cds)
        with self.assertRaisesRegex(ValueError, "not contained by"):
            self.cluster.add_cds(self.outside_cds)

    def test_mismatching_secmet_inside(self):
        self.add_core_function(self.inside_cds, cluster_product=False)
        self.cluster.add_cds(self.inside_cds)
        assert self.cluster.cds_children
        assert not self.cluster.definition_cdses

    def test_mismatching_secmet_neighbour(self):
        self.add_core_function(self.neighbour_cds, cluster_product=False)
        self.cluster.add_cds(self.neighbour_cds)
        assert self.cluster.cds_children
        assert not self.cluster.definition_cdses

    def test_mismatching_secmet_outside(self):
        self.add_core_function(self.outside_cds, cluster_product=False)
        with self.assertRaisesRegex(ValueError, "not contained by"):
            self.cluster.add_cds(self.outside_cds)


class TestConstruction(unittest.TestCase):
    def test_product(self):
        loc = FeatureLocation(1, 6, strand=1)
        for bad in ["-", "-like", "NRPS-", "NRPS PKS", "NRPS/PKS", "NRPS,PKS", "NRPS.PKS"]:
            with self.assertRaisesRegex(ValueError, "invalid protocluster product"):
                Protocluster(loc, loc, tool="test", cutoff=17, neighbourhood_range=5,
                             product=bad, detection_rule="some rule text")
