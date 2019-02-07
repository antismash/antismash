# License: GNU Affero General Public License v3 or later
# A copy of GNU AGPL v3 should have been included in this software package in LICENSE.txt.

# for test files, silence irrelevant and noisy pylint warnings
# pylint: disable=no-self-use,protected-access,missing-docstring

from tempfile import NamedTemporaryFile
import unittest

from Bio.Alphabet import generic_dna
from Bio.Seq import Seq
from helperlibs.bio import seqio

from antismash.common.secmet import FeatureLocation, Record
from antismash.common.secmet.features import CandidateCluster, SubRegion, Region
from antismash.common.secmet.test.helpers import (
    DummyCandidateCluster,
    DummyCDS,
    DummyProtocluster,
)


def create_protocluster(start, end, product='a'):
    return DummyProtocluster(start, end, product=product)


class TestRegionChildren(unittest.TestCase):
    def setUp(self):
        self.protocluster = DummyProtocluster()
        self.candidate = DummyCandidateCluster([self.protocluster])
        self.sub = SubRegion(self.protocluster.location, "testtool")
        self.region = Region(candidate_clusters=[self.candidate], subregions=[self.sub])

    def test_children_accessible(self):
        assert self.region.subregions == (self.sub,)
        assert self.region.candidate_clusters == (self.candidate,)

    def test_children_immutable(self):
        with self.assertRaisesRegex(AttributeError, "can't set attribute"):
            self.region.subregions = (self.candidate,)
        with self.assertRaisesRegex(AttributeError, "can't set attribute"):
            self.region.candidate_clusters = (self.sub,)
        with self.assertRaisesRegex(AttributeError, "can't set attribute"):
            self.region.cds_children = []

    def test_incorrect_args(self):
        with self.assertRaises(AssertionError):
            Region(candidate_clusters=[self.sub])
        with self.assertRaises(AssertionError):
            Region(subregions=[self.candidate])

    def test_missing_children(self):
        with self.assertRaisesRegex(ValueError, "at least one"):
            Region()
        with self.assertRaisesRegex(ValueError, "at least one"):
            Region(candidate_clusters=[], subregions=[])

    def test_add_cds_propagation(self):
        cds = DummyCDS(0, 10)
        assert cds.is_contained_by(self.region)
        # ensure all empty to start with
        assert not self.protocluster.cds_children
        assert not self.candidate.cds_children
        assert not self.sub.cds_children
        assert not self.region.cds_children
        assert not cds.region

        self.region.add_cds(cds)
        assert self.protocluster.cds_children == (cds,)
        assert self.candidate.cds_children == (cds,)
        assert self.sub.cds_children == (cds,)
        assert self.region.cds_children == (cds,)
        assert cds.region is self.region

    def test_limited_add_cds_propagation(self):
        cds = DummyCDS(0, 10)
        self.sub = SubRegion(FeatureLocation(20, 30), "testtool")
        self.region = Region(candidate_clusters=[self.candidate], subregions=[self.sub])

        # ensure all empty to start with
        assert not self.protocluster.cds_children
        assert not self.candidate.cds_children
        assert not self.sub.cds_children
        assert not self.region.cds_children
        assert not cds.region

        self.region.add_cds(cds)
        assert self.protocluster.cds_children == (cds,)
        assert self.candidate.cds_children == (cds,)
        assert not self.sub.cds_children
        assert self.region.cds_children == (cds,)
        assert cds.region is self.region

    def test_adding_invalid_cds(self):
        cds = DummyCDS(50, 60)
        assert not cds.is_contained_by(self.region)
        with self.assertRaisesRegex(ValueError, "not contained by"):
            self.region.add_cds(cds)

    def test_unique_clusters(self):
        protoclusters = [create_protocluster(i, 10, product=prod) for i, prod in enumerate("ABC")]
        candidates = [CandidateCluster(CandidateCluster.kinds.INTERLEAVED, protoclusters[:2]),
                      CandidateCluster(CandidateCluster.kinds.INTERLEAVED, protoclusters[1:])]
        assert protoclusters[1] in candidates[0].protoclusters and protoclusters[1] in candidates[1].protoclusters
        region = Region(candidate_clusters=candidates)
        unique_clusters = region.get_unique_protoclusters()
        # if the protocluster in both candidates is repeated, there'll be an extra
        assert len(unique_clusters) == 3
        assert unique_clusters == protoclusters


class TestRegion(unittest.TestCase):
    def test_products(self):
        candidates = [DummyCandidateCluster([create_protocluster(0, 10)])]
        region = Region(candidate_clusters=candidates)
        assert region.products == ["a"]
        assert region.get_product_string() == "a"

        candidates = []
        for i, prod in zip(range(2), "ba"):
            candidates.append(DummyCandidateCluster([create_protocluster(i*10, (i+1)*10, product=prod)]))
        region = Region(candidate_clusters=candidates)
        assert region.products == ["b", "a"]
        assert region.get_product_string() == "a,b"

    def test_probabilities(self):
        loc = FeatureLocation(0, 10)
        candidates = [DummyCandidateCluster([create_protocluster(0, 10)])]
        assert Region(candidate_clusters=candidates).probabilities == []
        subs = [SubRegion(loc, "testtool", probability=None)]
        assert Region(candidate_clusters=candidates, subregions=subs).probabilities == []
        subs.append(SubRegion(loc, "testtool", probability=0.1))
        assert Region(candidate_clusters=candidates, subregions=subs).probabilities == [0.1]
        subs.append(SubRegion(loc, "testtool", probability=0.7))
        assert Region(candidate_clusters=candidates, subregions=subs).probabilities == [0.1, 0.7]

    def test_genbank(self):
        dummy_record = Record(Seq("A"*100, generic_dna))
        clusters = [create_protocluster(3, 20, "prodA"),
                    create_protocluster(25, 41, "prodB")]
        for cluster in clusters:
            dummy_record.add_protocluster(cluster)
        subregion = SubRegion(FeatureLocation(35, 71), "test", 0.7)
        dummy_record.add_subregion(subregion)
        candidate = CandidateCluster(CandidateCluster.kinds.NEIGHBOURING, clusters)
        dummy_record.add_candidate_cluster(candidate)
        region = Region(candidate_clusters=[candidate],
                        subregions=[subregion])
        dummy_record.add_region(region)
        with NamedTemporaryFile(suffix=".gbk") as output:
            region.write_to_genbank(output.name)
            bio = list(seqio.parse(output.name))
        assert len(bio) == 1
        print(bio[0].features)
        rec = Record.from_biopython(bio[0], taxon="bacteria")
        assert len(rec.get_regions()) == 1
        new = rec.get_region(0)
        assert new.location.start == 3 - region.location.start
        assert new.location.end == 71 - region.location.start
        assert new.products == region.products
        assert new.probabilities == region.probabilities
