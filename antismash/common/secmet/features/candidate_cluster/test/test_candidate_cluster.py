# License: GNU Affero General Public License v3 or later
# A copy of GNU AGPL v3 should have been included in this software package in LICENSE.txt.

# for test files, silence irrelevant and noisy pylint warnings
# pylint: disable=use-implicit-booleaness-not-comparison,protected-access,missing-docstring

import unittest

from Bio.Seq import Seq

from antismash.common.secmet import FeatureLocation, Record
from antismash.common.secmet.features import Protocluster, CandidateCluster
from antismash.common.secmet.features.candidate_cluster import (
    CandidateClusterKind,
)
from antismash.common.secmet.qualifiers import GeneFunction
from antismash.common.secmet.test.helpers import DummyCDS


def create_cds(start, end, products):
    cds = DummyCDS(start, end, locus_tag=f"{start}-{end}-{'-'.join(products)}")
    for product in products:
        cds.gene_functions.add(GeneFunction.CORE, "test", "dummy", product)
    return cds


def create_cluster(n_start, start, end, n_end, product='a'):
    cluster = Protocluster(FeatureLocation(start, end),
                           FeatureLocation(n_start, n_end),
                           tool="testing", product=product, cutoff=1,
                           neighbourhood_range=0, detection_rule="some rule text")
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
        real = CandidateCluster.from_biopython(bio, record=self.record)
        assert isinstance(real, CandidateCluster)
        assert len(real.protoclusters) == len(self.record.get_protoclusters())
        for reference, record_cluster in zip(real.protoclusters, self.record.get_protoclusters()):
            assert reference is record_cluster

        # attempt a conversion with a record missing the cluster
        self.record.clear_protoclusters()
        with self.assertRaisesRegex(ValueError, "record does not contain all expected protoclusters"):
            CandidateCluster.from_biopython(bio, record=self.record)
        # and with no record
        with self.assertRaisesRegex(ValueError, "record instance required"):
            CandidateCluster.from_biopython(bio)

    def test_core(self):
        protos = [create_cluster(5, 10, 20, 25, "a"),
                  create_cluster(30, 40, 50, 60, "b")]
        cluster = CandidateCluster(CandidateClusterKind.NEIGHBOURING, protos,
                                   smiles="dummy", polymer="dummy")
        assert cluster.core_location == FeatureLocation(10, 50)

    def test_comparison(self):
        candidate = CandidateCluster(CandidateClusterKind.NEIGHBOURING, [create_cluster(5, 10, 20, 25, "a")])
        longer = CandidateCluster(CandidateClusterKind.NEIGHBOURING, [create_cluster(5, 10, 40, 45, "a")])
        after = CandidateCluster(CandidateClusterKind.NEIGHBOURING, [create_cluster(10, 20, 40, 45, "a")])

        def check(first, second):
            assert first < second
            assert first < second.location
            assert sorted([second, first]) == [first, second]

        check(candidate, after)
        check(longer, candidate)
        check(longer, after)
        assert sorted([after, candidate, longer]) == [longer, candidate, after]
