# License: GNU Affero General Public License v3 or later
# A copy of GNU AGPL v3 should have been included in this software package in LICENSE.txt.

# for test files, silence irrelevant and noisy pylint warnings
# pylint: disable=no-self-use,protected-access,missing-docstring

import unittest

from Bio.Seq import Seq

from antismash.common.secmet import FeatureLocation, Record
from antismash.common.secmet.features import Cluster, SuperCluster
from antismash.common.secmet.features.supercluster import (
    SuperClusterKind,
    TemporarySuperCluster,
    create_superclusters_from_clusters as creator,
)
from antismash.common.secmet.qualifiers import GeneFunction
from antismash.common.secmet.test.helpers import DummyCDS


def create_cds(start, end, products):
    cds = DummyCDS(start, end, locus_tag="%s-%s-%s" % (start, end, "-".join(products)))
    for product in products:
        cds.gene_functions.add(GeneFunction.CORE, "test", "dummy", product)
    return cds


def create_cluster(n_start, start, end, n_end, product='a'):
    cluster = Cluster(FeatureLocation(start, end),
                      FeatureLocation(n_start, n_end),
                      tool="testing", product=product, cutoff=1,
                      neighbourhood_range=0, detection_rule="some rule text")
    cds = create_cds(start, end, [product])
    cluster.add_cds(cds)
    return cluster


class TestKind(unittest.TestCase):
    def test_str_conversion(self):
        for kind in SuperClusterKind:
            assert SuperClusterKind.from_string(str(kind)) == kind
        with self.assertRaisesRegex(ValueError, "unknown"):
            SuperClusterKind.from_string("test")


class TestSuperCluster(unittest.TestCase):
    def setUp(self):
        self.record = Record(Seq("A"*100))
        clusters = [create_cluster(0, 0, 10, 10)]
        for cluster in clusters:
            self.record.add_cluster(cluster)

    def test_kinds_attachment(self):
        assert SuperCluster.kinds == SuperClusterKind

    def test_record_linkage(self):
        cluster = SuperCluster(SuperCluster.kinds.INTERLEAVED, self.record.get_clusters())
        with self.assertRaisesRegex(ValueError, "SuperCluster not contained in record"):
            cluster.get_supercluster_number()
        self.record.add_supercluster(cluster)
        assert cluster.get_supercluster_number() == 1

    def test_bad_kind(self):
        with self.assertRaisesRegex(TypeError, "should be SuperClusterKind"):
            SuperCluster("berf", self.record.get_clusters())

    def test_no_clusters(self):
        with self.assertRaisesRegex(ValueError, "cannot exist without at least one"):
            SuperCluster(SuperCluster.kinds.INTERLEAVED, [])

    def test_rules(self):
        cluster = SuperCluster(SuperCluster.kinds.INTERLEAVED, self.record.get_clusters())
        assert cluster.detection_rules == [cluster.detection_rule for cluster in self.record.get_clusters()]

    def test_smiles_and_polymer(self):
        cluster = SuperCluster(SuperCluster.kinds.INTERLEAVED, self.record.get_clusters())
        assert cluster.smiles_structure is None
        assert cluster.polymer is None

    def test_conversion(self):
        kind = SuperClusterKind.INTERLEAVED
        original = SuperCluster(kind, self.record.get_clusters(),
                                smiles="dummy smiles", polymer="dummy polymer")
        self.record.add_supercluster(original)
        assert original.products == ["a"]
        assert len(original.clusters) == 1
        bios = original.to_biopython()
        assert len(bios) == 1
        bio = bios[0]
        assert bio.qualifiers["product"] == ["a"]
        assert bio.qualifiers["kind"] == [str(kind)]
        assert bio.qualifiers["supercluster_number"] == [str(original.get_supercluster_number())]
        assert bio.qualifiers["SMILES"] == ["dummy smiles"]
        assert bio.qualifiers["polymer"] == ["dummy polymer"]
        assert bio.qualifiers["contig_edge"] == ["True"]
        regenerated = SuperCluster.from_biopython(bio)
        assert isinstance(regenerated, TemporarySuperCluster)
        assert regenerated.products == original.products
        assert regenerated.location == original.location
        assert regenerated.smiles_structure == original.smiles_structure
        assert regenerated.polymer == original.polymer
        assert regenerated.clusters == [cluster.get_cluster_number() for cluster in self.record.get_clusters()]
        assert regenerated.kind == original.kind

        real = regenerated.convert_to_real_feature(self.record)
        assert isinstance(real, SuperCluster)
        assert len(real.clusters) == len(self.record.get_clusters())
        for reference, record_cluster in zip(real.clusters, self.record.get_clusters()):
            assert reference is record_cluster

        # attempt a conversion with a record missing the cluster
        self.record.clear_clusters()
        with self.assertRaisesRegex(ValueError, "Not all referenced clusters are present"):
            regenerated.convert_to_real_feature(self.record)


class TestCreation(unittest.TestCase):
    def test_creation_empty(self):
        assert creator([]) == []

    def test_creation_single(self):
        created = creator([create_cluster(3, 8, 71, 76, 'a')])
        print(created)
        assert len(created) == 1
        assert created[0].location == FeatureLocation(3, 76)
        assert created[0].kind == SuperCluster.kinds.SINGLE

    def test_creation_neighbours(self):
        cluster = create_cluster(3, 8, 71, 76, 'a')
        extra_cluster = create_cluster(50, 100, 120, 170, 'b')
        created = creator([cluster, extra_cluster])
        print(created)
        assert len(created) == 3
        expected_location = FeatureLocation(cluster.location.start,
                                            extra_cluster.location.end)
        assert created[0].kind == SuperCluster.kinds.NEIGHBOURING and created[0].location == expected_location
        assert created[1].kind == SuperCluster.kinds.SINGLE and created[1].location == cluster.location
        assert created[2].kind == SuperCluster.kinds.SINGLE and created[2].location == extra_cluster.location

    def test_creation_coreoverlap(self):
        cluster = create_cluster(3, 8, 71, 76, 'a')
        extra_cluster = create_cluster(50, 60, 120, 170, 'b')
        # create a CDS within both clusters that has a product from only one cluster
        cds = create_cds(60, 65, ["a"])
        cluster.add_cds(cds)
        extra_cluster.add_cds(cds)

        created = creator([cluster, extra_cluster])
        print(created)
        assert len(created) == 1
        supercluster = created[0]
        assert supercluster.kind == SuperCluster.kinds.INTERLEAVED
        assert supercluster.location == FeatureLocation(3, 170)

    def test_creation_hybrid(self):
        cluster = create_cluster(3, 8, 71, 76, 'a')
        hybrid_cluster = create_cluster(50, 60, 120, 170, 'b')

        # insert the cds that will cause the hybrid call
        cds_ab = create_cds(60, 65, ["a", "b"])
        cluster.add_cds(cds_ab)
        hybrid_cluster.add_cds(cds_ab)

        created = creator([cluster, hybrid_cluster])
        print(created)
        assert len(created) == 1
        supercluster = created[0]
        assert supercluster.kind == SuperCluster.kinds.CHEMICAL_HYBRID
        assert supercluster.location == FeatureLocation(3, 170)

    def test_creation_mixed(self):
        cluster = create_cluster(3, 8, 71, 76, 'a')
        hybrid_cluster = create_cluster(50, 60, 120, 170, 'b')
        overlap_cluster = create_cluster(80, 90, 130, 180, 'o')
        neighbour_cluster = create_cluster(50, 210, 260, 270, 'a')
        isolated_cluster = create_cluster(450, 500, 550, 600, 'alone')

        # insert the cds that will cause the hybrid call
        cds_ab = create_cds(60, 65, ["a", "b"])
        cluster.add_cds(cds_ab)
        hybrid_cluster.add_cds(cds_ab)

        created = creator([cluster, hybrid_cluster, overlap_cluster, neighbour_cluster, isolated_cluster])
        print(created)

        assert len(created) == 5
        assert created[0].location == FeatureLocation(3, 270)
        assert created[0].kind == SuperCluster.kinds.NEIGHBOURING
        assert created[0].clusters == (cluster, hybrid_cluster, overlap_cluster, neighbour_cluster)

        assert created[1].location == FeatureLocation(3, 180)
        assert created[1].kind == SuperCluster.kinds.INTERLEAVED
        assert created[1].clusters == (cluster, hybrid_cluster, overlap_cluster)

        assert created[2].location == FeatureLocation(3, 170)
        assert created[2].kind == SuperCluster.kinds.CHEMICAL_HYBRID
        assert created[2].clusters == (cluster, hybrid_cluster)

        assert created[3].location == FeatureLocation(50, 270)
        assert created[3].kind == SuperCluster.kinds.SINGLE
        assert created[3].clusters == (neighbour_cluster,)

        assert created[4].location == FeatureLocation(450, 600)
        assert created[4].kind == SuperCluster.kinds.SINGLE
        assert created[4].clusters == (isolated_cluster,)
