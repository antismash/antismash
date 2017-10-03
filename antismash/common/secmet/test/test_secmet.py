# License: GNU Affero General Public License v3 or later
# A copy of GNU AGPL v3 should have been included in this software package in LICENSE.txt.

# for test files, silence irrelevant and noisy pylint warnings
# pylint: disable=no-self-use,protected-access,missing-docstring

from collections import defaultdict
import unittest

import Bio.SeqIO
from Bio.Seq import Seq
from Bio.SeqFeature import SeqFeature, FeatureLocation
from antismash.common.test import helpers
from antismash.common.secmet import Record, Cluster, CDSFeature, Feature, GeneFunction

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
            record.add_cluster(cluster)
        for i, cluster in enumerate(sorted(list(record.get_clusters()))):
            assert cluster.get_cluster_number() == i + 1

    def test_overlapping_clusters(self):
        record = Record()
        record.add_cluster(Cluster(FeatureLocation(10, 40), 0, 0, []))
        with self.assertRaises(ValueError):
            record.add_cluster(Cluster(FeatureLocation(0, 11), 0, 0, []))
        # ok, since ends aren't inclusive
        record.add_cluster(Cluster(FeatureLocation(0, 10), 0, 0, []))

    def test_cds_cluster_linkage(self):
        record = Record()
        for start, end in [(50, 100), (10, 90), (0, 9), (150, 200)]:
            record.add_cds_feature(helpers.DummyCDS(start, end))
        for start, end in [(10, 120), (5, 110), (10, 160), (45, 200)]:
            record.clear_clusters()
            cluster = helpers.DummyCluster(start, end)
            record.add_cluster(cluster)
            assert len(cluster.cds_children) == 2
            for cds in cluster.cds_children:
                assert cds.overlaps_with(cluster)

class TestFeature(unittest.TestCase):
    def test_overlaps_with(self):
        # no overlap
        a = helpers.DummyFeature(5, 10)
        assert isinstance(a, Feature) # just to be sure it works the way we want
        b = helpers.DummyFeature(100, 110)
        assert not a.overlaps_with(b) and not b.overlaps_with(a)
        # completely within
        b = helpers.DummyFeature(0, 20)
        assert a.overlaps_with(b) and b.overlaps_with(a)
        # partially within
        b = helpers.DummyFeature(0, 8)
        assert a.overlaps_with(b) and b.overlaps_with(a)
        # borders touching
        b = helpers.DummyFeature(0, 5)
        assert not (a.overlaps_with(b) or b.overlaps_with(a))

    def test_is_contained_by(self):
        # same location is considered to be contained
        a = helpers.DummyFeature(5, 10)
        assert a.is_contained_by(a)
        for strand in (-1, 1):
            # no overlap
            b = helpers.DummyFeature(15, 25, strand)
            assert not a.is_contained_by(b)
            assert not b.is_contained_by(a)
            # b is contained
            b = helpers.DummyFeature(6, 9, strand)
            assert not a.is_contained_by(b)
            assert b.is_contained_by(a)
            # only partial overlap
            b = helpers.DummyFeature(6, 19, strand)
            assert not a.is_contained_by(b)
            assert not b.is_contained_by(a)
            b = helpers.DummyFeature(1, 7, strand)
            assert not a.is_contained_by(b)
            assert not b.is_contained_by(a)
            # edge cases
            b = helpers.DummyFeature(5, 7, strand)
            assert not a.is_contained_by(b)
            assert b.is_contained_by(a)
            b = helpers.DummyFeature(7, 10, strand)
            assert not a.is_contained_by(b)
            assert b.is_contained_by(a)

    def test_biopython_conversion(self):
        bio = SeqFeature(FeatureLocation(1, 5))
        bio.qualifiers["foo"] = ["bar"]
        # check that features without types are caught
        with self.assertRaises(AssertionError):
            sec = Feature.from_biopython(bio)
        bio.type = "test"
        sec = Feature.from_biopython(bio)
        assert sec.get_qualifier("foo") == tuple(["bar"])
        assert sec.get_qualifier("bar") is None

    def test_string_conversion(self):
        for feature_type in ["cluster", "cds_motif", "test"]:
            for start, end in [(1, 5), (3, 8), (10, 15)]:
                feature = Feature(FeatureLocation(start, end, strand=1),
                                  feature_type=feature_type)
                assert str(feature) == "%s([%d:%d](+))" % (feature_type, start, end)
                feature = Feature(FeatureLocation(start, end, strand=-1),
                                  feature_type=feature_type)
                assert str(feature) == "%s([%d:%d](-))" % (feature_type, start, end)

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

class TestCDSFeature(unittest.TestCase):
    def test_required_identifiers(self):
        with self.assertRaises(ValueError):
            dummy = CDSFeature(FeatureLocation(1, 5))
        dummy = CDSFeature(FeatureLocation(1, 5), locus_tag="foo")
        dummy = CDSFeature(FeatureLocation(1, 5), protein_id="foo")
        dummy = CDSFeature(FeatureLocation(1, 5), gene="foo")

class TestGeneFunction(unittest.TestCase):
    def test_membership(self):
        dummy = GeneFunction.OTHER
        with self.assertRaises(AttributeError):
            dummy = GeneFunction.non_existant

    def test_equality(self):
        assert GeneFunction.OTHER == GeneFunction.OTHER
        assert GeneFunction.CORE != GeneFunction.OTHER

    def test_string_conversion(self):
        assert str(GeneFunction.OTHER) == "other"
        for member in dir(GeneFunction):
            if member.isupper():
                assert str(getattr(GeneFunction, member)) == member.lower()

    def test_CDS_function(self):
        cds = CDSFeature(FeatureLocation(1, 5), locus_tag="foo")
        # default value
        assert cds.gene_function == GeneFunction.OTHER
        # check bad values can't be assigned
        with self.assertRaises(AssertionError):
            cds.gene_function = "other"
        with self.assertRaises(AssertionError):
            cds.gene_function = 0
        # check overriding OTHER works
        cds.gene_function = GeneFunction.CORE
        assert cds.gene_function == GeneFunction.CORE
        # test that overriding non-other doesn't work
        cds.gene_function = GeneFunction.ADDITIONAL
        assert cds.gene_function == GeneFunction.CORE
