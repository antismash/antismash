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
from antismash.common.secmet.features import Cluster, SuperCluster, SubRegion, Region
from antismash.common.secmet.test.helpers import DummyCDS


def create_cluster(start, end, product='a'):
    return Cluster(FeatureLocation(start, end),
                   FeatureLocation(start, end),
                   tool="testing", product=product, cutoff=1,
                   neighbourhood_range=0, detection_rule="some rule text")


class TestRegionChildren(unittest.TestCase):
    def setUp(self):
        self.cluster = create_cluster(0, 10)
        self.super = SuperCluster(SuperCluster.kinds.SINGLE, [self.cluster])
        self.sub = SubRegion(self.cluster.location, "testtool")
        self.region = Region(superclusters=[self.super], subregions=[self.sub])

    def test_children_accessible(self):
        assert self.region.subregions == (self.sub,)
        assert self.region.superclusters == (self.super,)

    def test_children_immutable(self):
        with self.assertRaisesRegex(AttributeError, "can't set attribute"):
            self.region.subregions = (self.super,)
        with self.assertRaisesRegex(AttributeError, "can't set attribute"):
            self.region.superclusters = (self.sub,)
        with self.assertRaisesRegex(AttributeError, "can't set attribute"):
            self.region.cds_children = []

    def test_incorrect_args(self):
        with self.assertRaises(AssertionError):
            Region(superclusters=[self.sub])
        with self.assertRaises(AssertionError):
            Region(subregions=[self.super])

    def test_missing_children(self):
        with self.assertRaisesRegex(ValueError, "at least one"):
            Region()
        with self.assertRaisesRegex(ValueError, "at least one"):
            Region(superclusters=[], subregions=[])

    def test_add_cds_propagation(self):
        cds = DummyCDS(0, 10)
        assert cds.is_contained_by(self.region)
        # ensure all empty to start with
        assert not self.cluster.cds_children
        assert not self.super.cds_children
        assert not self.sub.cds_children
        assert not self.region.cds_children
        assert not cds.region

        self.region.add_cds(cds)
        assert self.cluster.cds_children == (cds,)
        assert self.super.cds_children == (cds,)
        assert self.sub.cds_children == (cds,)
        assert self.region.cds_children == (cds,)
        assert cds.region is self.region

    def test_limited_add_cds_propagation(self):
        cds = DummyCDS(0, 10)
        self.sub = SubRegion(FeatureLocation(20, 30), "testtool")
        self.region = Region(superclusters=[self.super], subregions=[self.sub])

        # ensure all empty to start with
        assert not self.cluster.cds_children
        assert not self.super.cds_children
        assert not self.sub.cds_children
        assert not self.region.cds_children
        assert not cds.region

        self.region.add_cds(cds)
        assert self.cluster.cds_children == (cds,)
        assert self.super.cds_children == (cds,)
        assert not self.sub.cds_children
        assert self.region.cds_children == (cds,)
        assert cds.region is self.region

    def test_adding_invalid_cds(self):
        cds = DummyCDS(50, 60)
        assert not cds.is_contained_by(self.region)
        with self.assertRaisesRegex(ValueError, "not contained by"):
            self.region.add_cds(cds)

    def test_unique_clusters(self):
        clusters = [create_cluster(i, 10, product=prod) for i, prod in enumerate("ABC")]
        superclusters = [SuperCluster(SuperCluster.kinds.INTERLEAVED, clusters[:2]),
                         SuperCluster(SuperCluster.kinds.INTERLEAVED, clusters[1:])]
        assert clusters[1] in superclusters[0].clusters and clusters[1] in superclusters[1].clusters
        region = Region(superclusters=superclusters)
        unique_clusters = region.get_unique_clusters()
        # if the cluster in both superclusters is repeated, there'll be an extra
        assert len(unique_clusters) == 3
        assert unique_clusters == clusters


class TestRegion(unittest.TestCase):
    def test_products(self):
        supers = [SuperCluster(SuperCluster.kinds.SINGLE, [create_cluster(0, 10)])]
        region = Region(superclusters=supers)
        assert region.products == ["a"]
        assert region.get_product_string() == "a"

        supers = []
        for i, prod in zip(range(2), "ba"):
            supers.append(SuperCluster(SuperCluster.kinds.SINGLE, [create_cluster(i*10, (i+1)*10, product=prod)]))
        region = Region(superclusters=supers)
        assert region.products == ["b", "a"]
        assert region.get_product_string() == "a-b"

    def test_probabilities(self):
        loc = FeatureLocation(0, 10)
        supers = [SuperCluster(SuperCluster.kinds.SINGLE, [create_cluster(0, 10)])]
        assert Region(superclusters=supers).probabilities == []
        subs = [SubRegion(loc, "testtool", probability=None)]
        assert Region(superclusters=supers, subregions=subs).probabilities == []
        subs.append(SubRegion(loc, "testtool", probability=0.1))
        assert Region(superclusters=supers, subregions=subs).probabilities == [0.1]
        subs.append(SubRegion(loc, "testtool", probability=0.7))
        assert Region(superclusters=supers, subregions=subs).probabilities == [0.1, 0.7]

    def test_genbank(self):
        dummy_record = Record(Seq("A"*100, generic_dna))
        clusters = [create_cluster(3, 20, "prodA"),
                    create_cluster(25, 41, "prodB")]
        for cluster in clusters:
            dummy_record.add_cluster(cluster)
        subregion = SubRegion(FeatureLocation(35, 71), "test", 0.7)
        dummy_record.add_subregion(subregion)
        supercluster = SuperCluster(SuperCluster.kinds.NEIGHBOURING, clusters)
        dummy_record.add_supercluster(supercluster)
        region = Region(superclusters=[supercluster],
                        subregions=[subregion])
        dummy_record.add_region(region)
        with NamedTemporaryFile(suffix=".gbk") as output:
            region.write_to_genbank(output.name)
            bio = list(seqio.parse(output.name))
        assert len(bio) == 1
        rec = Record.from_biopython(bio[0], taxon="bacteria")
        assert len(rec.get_regions()) == 1
        new = rec.get_region(0)
        assert new.location.start == 3 - region.location.start
        assert new.location.end == 71 - region.location.start
        assert new.products == region.products
        assert new.probabilities == region.probabilities
