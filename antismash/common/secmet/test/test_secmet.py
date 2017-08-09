# License: GNU Affero General Public License v3 or later
# A copy of GNU AGPL v3 should have been included in this software package in LICENSE.txt.

from collections import defaultdict
import unittest

import Bio.SeqIO
from Bio.Seq import Seq
from Bio.SeqFeature import FeatureLocation
from antismash.common.test import helpers
from antismash.common.secmet import Record, Cluster, CDSFeature, Feature

class TestConversion(unittest.TestCase):
    def test_conversion(self):
        before = list(Bio.SeqIO.parse(helpers.get_path_to_nisin_genbank(), "genbank"))[0]
        # sort notes, because direct comparisons otherwise are awful
        for feature in before.features:
            if "note" in feature.qualifiers:
                feature.qualifiers["note"] = sorted(feature.qualifiers["note"])
        before_features = sorted(map(str, before.features))
        type_counts = defaultdict(lambda: 0)
        for feature in before.features:
            type_counts[feature.type] += 1
        record = Record.from_biopython(before)
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
        for i, cluster in enumerate(sorted(list(record.get_clusters()))):
            assert cluster.get_cluster_number() == i

    def test_overlapping_clusters(self):
        record = Record()
        record.add_cluster(Cluster(FeatureLocation(10, 40), 0, 0, []))
        with self.assertRaises(ValueError):
            record.add_cluster(Cluster(FeatureLocation(0, 10), 0, 0, []))

    def test_cds_cluster_linkage(self):
        record = Record()
        for start, end in [(50, 100), (10, 90), (0, 9), (150, 200)]:
            record.add_cds_feature(helpers.DummyCDS(start, end))
        for start, end in [(25, 120), (10, 90), (10, 140), (80, 140), (95, 190)]:
            record.clear_clusters()
            cluster = helpers.DummyCluster(start, end)
            record.add_cluster(cluster)
            assert len(cluster.cds_children) == 2
            for cds in cluster.cds_children:
                assert cds.overlaps_with(cluster)

class TestFeature(unittest.TestCase):
    def test_overlaps_with(self):
        # no overlap
        a = Feature(FeatureLocation(5, 10), feature_type="test")
        b = Feature(FeatureLocation(100, 110), feature_type="test")
        assert not a.overlaps_with(b) and not b.overlaps_with(a)
        # completely within
        b = Feature(FeatureLocation(0, 20), feature_type="test")
        assert a.overlaps_with(b) and b.overlaps_with(a)
        # partially within
        b = Feature(FeatureLocation(0, 8), feature_type="test")
        assert a.overlaps_with(b) and b.overlaps_with(a)
        # borders touching
        b = Feature(FeatureLocation(0, 5), feature_type="test")
        assert a.overlaps_with(b) and b.overlaps_with(a)

# since we're about to test assigning to non-slots, shut pylint up
# pylint: disable=assigning-non-slot
    def test_membership(self):
        location = FeatureLocation(0, 1)
        # Features don't have locus tags
        with self.assertRaises(AttributeError):
            Feature(location, feature_type="none").locus_tag = "something"
        # CDSFeatures don't have an 'other_value'
        with self.assertRaises(AttributeError):
            CDSFeature(location, translation="none", gene="a").other_value = 1
        cluster = Cluster(location, 0, 0, products=["a", "b"])
        assert cluster.products == ["a", "b"]
        # Clusters have products, not product
        with self.assertRaises(AttributeError):
            cluster.product = ["c", "d"]
# pylint: enable=assigning-non-slot
