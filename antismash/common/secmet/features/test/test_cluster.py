# License: GNU Affero General Public License v3 or later
# A copy of GNU AGPL v3 should have been included in this software package in LICENSE.txt.

# for test files, silence irrelevant and noisy pylint warnings
# pylint: disable=no-self-use,protected-access,missing-docstring

from tempfile import NamedTemporaryFile
import unittest

from Bio.Alphabet import generic_dna
from Bio.Seq import Seq
from helperlibs.bio import seqio

from antismash.common.secmet import CDSFeature, FeatureLocation, Record
from antismash.common.secmet.features.cluster import Cluster


class TestCluster(unittest.TestCase):
    def create_cluster(self, start, end):
        return Cluster(FeatureLocation(start, end, strand=1),
                       cutoff=1, extent=1, products=['a'])

    def create_cds(self, start, end, strand=1):
        return CDSFeature(FeatureLocation(start, end, strand),
                          locus_tag="%d-%d" % (start, end))

    def setUp(self):
        self.record = Record(Seq("A" * 1000))
        self.start = 100
        self.end = 900
        self.cluster = self.create_cluster(self.start, self.end)
        self.record.add_cluster(self.cluster)
        assert self.cluster.location.start == self.start
        assert self.cluster.location.end == self.end

    def test_trim_unattached(self):
        cluster = self.create_cluster(1, 2)
        cluster.trim_overlapping()
        assert self.cluster.location.start == self.start
        assert self.cluster.location.end == self.end

    def test_trim_empty(self):
        self.cluster.trim_overlapping()
        assert self.cluster.location.start == self.start
        assert self.cluster.location.end == self.end

    def test_trim_contained(self):
        starts = [200, 300, 500]
        ends = [250, 350, 600]
        for start, end in zip(starts, ends):
            feature = self.create_cds(start, end)
            self.record.add_cds_feature(feature)
            assert feature.cluster == self.cluster
        self.cluster.trim_overlapping()
        assert self.cluster.location.start == self.start
        assert self.cluster.location.end == self.end

        for cds in self.record.get_cds_features():
            self.record.remove_cds_feature(cds)
        assert not self.cluster.cds_children

        for start, end in zip(starts, ends):
            feature = self.create_cds(start, end, strand=-1)
            self.record.add_cds_feature(feature)
            assert feature.cluster == self.cluster
        self.cluster.trim_overlapping()
        assert self.cluster.location.start == self.start
        assert self.cluster.location.end == self.end

    def test_trim_leading_overlap(self):
        self.record.add_cds_feature(self.create_cds(self.start - 3, self.start + 3))
        self.record.add_cds_feature(self.create_cds(self.start + 20, self.end - 20))
        self.cluster.trim_overlapping()
        assert self.cluster.location.start == self.start + 3
        assert self.cluster.location.end == self.end

    def test_trim_leading_overlap_with_overlapping_contained(self):  # pylint: disable=invalid-name
        self.record.add_cds_feature(self.create_cds(self.start - 3, self.start + 3))
        self.record.add_cds_feature(self.create_cds(self.start + 1, self.start + 10))
        self.cluster.trim_overlapping()
        assert self.cluster.location.start == self.start + 1
        assert self.cluster.location.end == self.end

    def test_trim_trailing_overlap(self):
        self.record.add_cds_feature(self.create_cds(self.end - 3, self.end + 3))
        self.record.add_cds_feature(self.create_cds(self.start + 20, self.end - 20))
        self.cluster.trim_overlapping()
        assert self.cluster.location.start == self.start
        assert self.cluster.location.end == self.end - 3

    def test_trim_trailing_overlap_with_overlapping_contained(self):  # pylint: disable=invalid-name
        self.record.add_cds_feature(self.create_cds(self.end - 3, self.end + 3))
        self.record.add_cds_feature(self.create_cds(self.end - 10, self.end - 1))
        self.cluster.trim_overlapping()
        assert self.cluster.location.start == self.start
        assert self.cluster.location.end == self.end - 1

    def test_products(self):
        assert self.cluster.products == ("a",)
        self.cluster.add_product("b")
        assert self.cluster.products == ("a", "b")
        with self.assertRaises(AttributeError):
            self.cluster.products.append("c")  # pylint: disable=no-member
        with self.assertRaises(AssertionError):
            self.cluster.add_product(None)
        with self.assertRaises(AssertionError):
            self.cluster.add_product(["C"])


class TestConversion(unittest.TestCase):
    def test_core(self):
        cluster = Cluster(FeatureLocation(3, 71, strand=1),
                          cutoff=17, extent=5, products=['a', 'c'])

        bio = cluster.to_biopython()
        assert len(bio) == 1
        new = Cluster.from_biopython(bio[0])
        assert new is not cluster
        assert new.cutoff == cluster.cutoff == 17
        assert new.extent == cluster.extent == 5
        assert new.products == cluster.products == ('a', 'c')
        assert new.location.start == cluster.location.start == 3

    def test_genbank(self):
        dummy_record = Record(Seq("A"*100, generic_dna))
        cluster = Cluster(FeatureLocation(3, 71, strand=1),
                          cutoff=17, extent=5, products=['a', 'c'])
        dummy_record.add_cluster(cluster)
        with NamedTemporaryFile(suffix=".gbk") as output:
            cluster.write_to_genbank(output.name)
            bio = list(seqio.parse(output.name))
        assert len(bio) == 1
        rec = Record.from_biopython(bio[0], taxon="bacteria")
        assert len(rec.get_clusters()) == 1
        new = rec.get_cluster(0)
        assert new.cutoff == cluster.cutoff
        assert new.products == cluster.products
        assert new.location.start == 0
        assert new.location.end == cluster.location.end - cluster.location.start
