# License: GNU Affero General Public License v3 or later
# A copy of GNU AGPL v3 should have been included in this software package in LICENSE.txt.

# for test files, silence irrelevant and noisy pylint warnings
# pylint: disable=no-self-use,protected-access,missing-docstring

import unittest

from Bio.Seq import Seq

from antismash.common.secmet import FeatureLocation, Record
from antismash.common.secmet.features import Protocluster, CandidateCluster
from antismash.common.secmet.features.candidate_cluster import (
    CandidateClusterKind,
    TemporaryCandidateCluster,
    create_candidates_from_protoclusters as creator,
)
from antismash.common.secmet.qualifiers import GeneFunction
from antismash.common.secmet.test.helpers import DummyCDS


def create_cds(start, end, products):
    cds = DummyCDS(start, end, locus_tag="%s-%s-%s" % (start, end, "-".join(products)))
    for product in products:
        cds.gene_functions.add(GeneFunction.CORE, "test", "dummy", product)
    return cds


def create_cluster(n_start, start, end, n_end, product='a'):
    cluster = Protocluster(FeatureLocation(start, end),
                           FeatureLocation(n_start, n_end),
                           tool="testing", product=product, cutoff=1,
                           neighbourhood_range=0, detection_rule="some rule text")
    cds = create_cds(start, end, [product])
    cluster.add_cds(cds)
    return cluster


class TestKind(unittest.TestCase):
    def test_str_conversion(self):
        for kind in CandidateClusterKind:
            assert CandidateClusterKind.from_string(str(kind)) == kind
        with self.assertRaisesRegex(ValueError, "unknown"):
            CandidateClusterKind.from_string("test")


class TestCandidateCluster(unittest.TestCase):
    def setUp(self):
        self.record = Record(Seq("A"*100))
        clusters = [create_cluster(0, 0, 10, 10)]
        for cluster in clusters:
            self.record.add_protocluster(cluster)

    def test_kinds_attachment(self):
        assert CandidateCluster.kinds == CandidateClusterKind

    def test_record_linkage(self):
        cluster = CandidateCluster(CandidateCluster.kinds.INTERLEAVED, self.record.get_protoclusters())
        with self.assertRaisesRegex(ValueError, "CandidateCluster not contained in record"):
            cluster.get_candidate_cluster_number()
        self.record.add_candidate_cluster(cluster)
        assert cluster.get_candidate_cluster_number() == 1

    def test_bad_kind(self):
        with self.assertRaisesRegex(TypeError, "should be CandidateClusterKind"):
            CandidateCluster("berf", self.record.get_protoclusters())

    def test_no_clusters(self):
        with self.assertRaisesRegex(ValueError, "cannot exist without at least one"):
            CandidateCluster(CandidateCluster.kinds.INTERLEAVED, [])

    def test_rules(self):
        cluster = CandidateCluster(CandidateCluster.kinds.INTERLEAVED, self.record.get_protoclusters())
        assert cluster.detection_rules == [cluster.detection_rule for cluster in self.record.get_protoclusters()]

    def test_smiles_and_polymer(self):
        cluster = CandidateCluster(CandidateCluster.kinds.INTERLEAVED, self.record.get_protoclusters())
        assert cluster.smiles_structure is None
        assert cluster.polymer is None

    def test_conversion(self):
        kind = CandidateClusterKind.INTERLEAVED
        original = CandidateCluster(kind, self.record.get_protoclusters(),
                                    smiles="dummy smiles", polymer="dummy polymer")
        self.record.add_candidate_cluster(original)
        assert original.products == ["a"]
        assert len(original.protoclusters) == 1
        bios = original.to_biopython()
        assert len(bios) == 1
        bio = bios[0]
        assert bio.qualifiers["product"] == ["a"]
        assert bio.qualifiers["kind"] == [str(kind)]
        assert bio.qualifiers["candidate_cluster_number"] == [str(original.get_candidate_cluster_number())]
        assert bio.qualifiers["SMILES"] == ["dummy smiles"]
        assert bio.qualifiers["polymer"] == ["dummy polymer"]
        assert bio.qualifiers["contig_edge"] == ["True"]
        regenerated = CandidateCluster.from_biopython(bio)
        assert isinstance(regenerated, TemporaryCandidateCluster)
        assert regenerated.products == original.products
        assert regenerated.location == original.location
        assert regenerated.smiles_structure == original.smiles_structure
        assert regenerated.polymer == original.polymer
        proto_numbers = [cluster.get_protocluster_number() for cluster in self.record.get_protoclusters()]
        assert regenerated.protoclusters == proto_numbers
        assert regenerated.kind == original.kind

        real = regenerated.convert_to_real_feature(self.record)
        assert isinstance(real, CandidateCluster)
        assert len(real.protoclusters) == len(self.record.get_protoclusters())
        for reference, record_cluster in zip(real.protoclusters, self.record.get_protoclusters()):
            assert reference is record_cluster

        # attempt a conversion with a record missing the cluster
        self.record.clear_protoclusters()
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
        assert created[0].kind == CandidateCluster.kinds.SINGLE

    def test_creation_neighbours(self):
        cluster = create_cluster(3, 8, 71, 76, 'a')
        extra_cluster = create_cluster(50, 100, 120, 170, 'b')
        created = creator([cluster, extra_cluster])
        print(created)
        assert len(created) == 3
        expected_location = FeatureLocation(cluster.location.start,
                                            extra_cluster.location.end)
        assert created[0].kind == CandidateCluster.kinds.NEIGHBOURING and created[0].location == expected_location
        assert created[1].kind == CandidateCluster.kinds.SINGLE and created[1].location == cluster.location
        assert created[2].kind == CandidateCluster.kinds.SINGLE and created[2].location == extra_cluster.location

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
        candidate_cluster = created[0]
        assert candidate_cluster.kind == CandidateCluster.kinds.INTERLEAVED
        assert candidate_cluster.location == FeatureLocation(3, 170)

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
        candidate_cluster = created[0]
        assert candidate_cluster.kind == CandidateCluster.kinds.CHEMICAL_HYBRID
        assert candidate_cluster.location == FeatureLocation(3, 170)

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
        assert created[0].kind == CandidateCluster.kinds.NEIGHBOURING
        assert created[0].protoclusters == (cluster, hybrid_cluster, overlap_cluster, neighbour_cluster)

        assert created[1].location == FeatureLocation(3, 180)
        assert created[1].kind == CandidateCluster.kinds.INTERLEAVED
        assert created[1].protoclusters == (cluster, hybrid_cluster, overlap_cluster)

        assert created[2].location == FeatureLocation(3, 170)
        assert created[2].kind == CandidateCluster.kinds.CHEMICAL_HYBRID
        assert created[2].protoclusters == (cluster, hybrid_cluster)

        assert created[3].location == FeatureLocation(50, 270)
        assert created[3].kind == CandidateCluster.kinds.SINGLE
        assert created[3].protoclusters == (neighbour_cluster,)

        assert created[4].location == FeatureLocation(450, 600)
        assert created[4].kind == CandidateCluster.kinds.SINGLE
        assert created[4].protoclusters == (isolated_cluster,)
