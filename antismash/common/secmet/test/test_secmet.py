# License: GNU Affero General Public License v3 or later
# A copy of GNU AGPL v3 should have been included in this software package in LICENSE.txt.

# for test files, silence irrelevant and noisy pylint warnings
# pylint: disable=no-self-use,protected-access,missing-docstring

from collections import defaultdict
import unittest

import Bio.SeqIO
from Bio.Seq import Seq
from Bio.SeqFeature import SeqFeature, FeatureLocation, CompoundLocation
from antismash.common.test import helpers
from antismash.common import secmet
from antismash.common.secmet.feature import Cluster, CDSFeature, Feature, ClusterBorder, PFAMDomain
from antismash.common.secmet.qualifiers import GeneFunction, GOQualifier
from antismash.common.secmet.record import _build_products_from_borders, Record


class TestConversion(unittest.TestCase):
    def test_record_conversion_from_biopython(self):
        before = list(Bio.SeqIO.parse(helpers.get_path_to_nisin_genbank(), "genbank"))[0]
        # sort notes, because direct comparisons otherwise are awful
        for feature in before.features:
            if "note" in feature.qualifiers:
                feature.qualifiers["note"] = sorted(feature.qualifiers["note"])
        before_features = sorted(map(str, before.features))
        type_counts = defaultdict(lambda: 0)
        for feature in before.features:
            type_counts[feature.type] += 1
        record = Record.from_biopython(before, taxon="bacteria")
        after = record.to_biopython()

        # ensure new features are correct
        assert len(before_features) == len(after.features)
        for bef, aft in zip(before_features, sorted(map(str, after.features))):
            assert bef == aft

        # ensure we haven't changed the original record or feature list
        assert id(before) != id(after)
        assert id(before.features) != id(after.features)
        for i in range(len(before.features)):
            assert id(before.features[i]) != id(after.features[i])
        for bef, aft in zip(before_features, sorted(map(str, before.features))):
            assert bef == aft

        # ensure that the counts of each match
        assert type_counts["CDS"] == len(record.get_cds_features())
        assert type_counts["PFAM_domain"] == len(record.get_pfam_domains())
        assert type_counts["cluster"] == len(record.get_clusters())
        assert type_counts["aSDomain"] == len(record.get_antismash_domains())


class TestRecord(unittest.TestCase):
    def test_cluster_numbering(self):
        record = Record(Seq("A"*1000))
        for start, end in [(50, 100), (10, 40), (700, 1000), (0, 9)]:
            cluster = helpers.DummyCluster(start, end)
            record.add_cluster(cluster)
        for i, cluster in enumerate(sorted(list(record.get_clusters()))):
            assert cluster.get_cluster_number() == i + 1

    def test_overlapping_clusters(self):
        record = Record(seq="A"*40)
        record.add_cluster(Cluster(FeatureLocation(10, 40), 0, 0, []))
        with self.assertRaises(ValueError):
            record.add_cluster(Cluster(FeatureLocation(0, 11), 0, 0, []))
        # ok, since ends aren't inclusive
        record.add_cluster(Cluster(FeatureLocation(0, 10), 0, 0, []))

    def test_cds_cluster_linkage(self):
        record = Record("A"*200)
        for start, end in [(50, 100), (10, 90), (0, 9), (150, 200)]:
            record.add_cds_feature(helpers.DummyCDS(start, end))
        for start, end in [(10, 120), (5, 110), (10, 160), (45, 200)]:
            record.clear_clusters()
            cluster = helpers.DummyCluster(start, end)
            record.add_cluster(cluster)
            assert len(cluster.cds_children) == 2
            for cds in cluster.cds_children:
                assert cds.overlaps_with(cluster)

    def test_cds_removal(self):
        record = Record(Seq("A" * 1000))
        cluster = helpers.DummyCluster(0, 1000)
        record.add_cluster(cluster)

        first_cds = helpers.DummyCDS(0, 100, locus_tag="A")
        second_cds = helpers.DummyCDS(200, 300, locus_tag="B")
        record.add_cds_feature(first_cds)
        record.add_cds_feature(second_cds)

        assert len(record.get_cds_features()) == 2
        assert len(cluster.cds_children) == 2

        record.remove_cds_feature(first_cds)

        assert len(record.get_cds_features()) == 1
        assert len(cluster.cds_children) == 1
        assert record.get_cds_features()[0] is list(cluster.cds_children)[0]
        assert record.get_cds_features()[0].locus_tag == "B"

    def test_orphaned_cluster_number(self):
        record = Record(Seq("A" * 1000))
        cluster = helpers.DummyCluster(0, 1000)

        with self.assertRaisesRegex(ValueError, "Cluster not contained in record"):
            print(record.get_cluster_number(cluster))

        with self.assertRaisesRegex(ValueError, "Cluster not contained in record"):
            print(cluster.get_cluster_number())

    def test_products_from_borders(self):
        location = FeatureLocation(1, 10)
        border1 = ClusterBorder(location, "toolA", product="A")
        assert border1.high_priority_product
        border2 = ClusterBorder(location, "toolB", product="B")
        assert border1.high_priority_product
        assert _build_products_from_borders([border1, border2]) == ["A", "B"]
        assert _build_products_from_borders([border2, border1]) == ["B", "A"]
        border2 = ClusterBorder(location, "toolB", product="B", high_priority_product=False)
        assert _build_products_from_borders([border1, border2]) == ["A"]
        assert _build_products_from_borders([border2, border1]) == ["A"]
        border1.high_priority_product = False
        assert _build_products_from_borders([border1, border2]) == ["A", "B"]
        assert _build_products_from_borders([border2, border1]) == ["B", "A"]


# since we're about to test assigning to non-slots, shut pylint up
# pylint: disable=assigning-non-slot
    def test_membership(self):
        location = FeatureLocation(0, 1, strand=1)
        # Features don't have locus tags
        with self.assertRaises(AttributeError):
            Feature(location, feature_type="none").locus_tag = "something"
        # CDSFeatures don't have an 'other_value'
        with self.assertRaises(AttributeError):
            CDSFeature(location, translation="none", gene="a").other_value = 1
        cluster = Cluster(location, 0, 0, products=["a", "b"])
        assert cluster.products == ("a", "b")
        # Clusters have products, not product
        with self.assertRaises(AttributeError):
            cluster.product = ["c", "d"]
# pylint: enable=assigning-non-slot
