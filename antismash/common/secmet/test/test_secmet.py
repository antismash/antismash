# License: GNU Affero General Public License v3 or later
# A copy of GNU AGPL v3 should have been included in this software package in LICENSE.txt.

# for test files, silence irrelevant and noisy pylint warnings
# pylint: disable=no-self-use,protected-access,missing-docstring

from collections import defaultdict
import unittest

import Bio.SeqIO
from Bio.Alphabet.IUPAC import IUPACProtein
from Bio.Seq import Seq
from Bio.SeqFeature import FeatureLocation, SeqFeature

from antismash.common.test import helpers
from ..features import (
    CDSFeature,
    Cluster,
    Feature,
    Region,
    SubRegion,
    SuperCluster,
)
from ..errors import SecmetInvalidInputError
from .helpers import DummyCDS
from ..record import Record


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

    def test_protein_sequences_caught(self):
        before = list(Bio.SeqIO.parse(helpers.get_path_to_nisin_genbank(), "genbank"))[0]

        # as a sanity check, make sure it's a seq and it functions as expected
        assert isinstance(before.seq, Seq)
        Record.from_biopython(before, taxon="bacteria")

        before.seq = Seq("AAAA", IUPACProtein())
        with self.assertRaisesRegex(ValueError, "protein records are not supported"):
            Record.from_biopython(before, taxon="bacteria")

    def test_missing_locations_caught(self):
        rec = list(Bio.SeqIO.parse(helpers.get_path_to_nisin_genbank(), "genbank"))[0]
        Record.from_biopython(rec, taxon="bacteria")
        rec.features.append(SeqFeature(None, type="broken"))
        with self.assertRaisesRegex(SecmetInvalidInputError, "feature is missing location"):
            Record.from_biopython(rec, taxon="bacteria")


class TestRecordFeatureNumbering(unittest.TestCase):
    def setUp(self):
        self.pairs = [(50, 100), (10, 40), (700, 1000), (0, 9)]
        self.locations = [FeatureLocation(start, end) for start, end in self.pairs]
        self.record = Record(Seq("A"*1000))

    def test_cluster_numbering(self):
        features = []
        for start, end in self.pairs:
            cluster = helpers.DummyCluster(start, end)
            self.record.add_cluster(cluster)
            features.append(cluster)
        features = sorted(features)
        for i, cluster in enumerate(self.record.get_clusters()):
            assert cluster.get_cluster_number() == i + 1
            assert self.record.get_cluster(i + 1) is features[i]

    def test_supercluster_numbering(self):
        features = []
        for location in self.locations:
            supercluster = SuperCluster(SuperCluster.kinds.SINGLE,
                                        [helpers.DummyCluster(location.start, location.end)])
            self.record.add_supercluster(supercluster)
            features.append(supercluster)
        features = sorted(features)
        for i, cluster in enumerate(self.record.get_superclusters()):
            assert cluster.get_supercluster_number() == i + 1
            assert self.record.get_supercluster(i + 1) is features[i]

    def test_subregion_numbering(self):
        features = []
        for location in self.locations:
            region = SubRegion(location, "test")
            self.record.add_subregion(region)
            features.append(region)
        features = sorted(features)
        for i, cluster in enumerate(self.record.get_subregions()):
            assert cluster.get_subregion_number() == i + 1
            assert self.record.get_subregion(i + 1) is features[i]

    def test_region_numbering(self):
        features = []
        for location in self.locations:
            region = Region(subregions=[SubRegion(location, "test")])
            self.record.add_region(region)
            features.append(region)
        features = sorted(features)
        for i, cluster in enumerate(self.record.get_regions()):
            assert cluster.get_region_number() == i + 1
            assert self.record.get_region(i + 1) is features[i]


class TestRecord(unittest.TestCase):
    def test_overlapping_regions(self):
        record = Record("A"*40)
        record.add_region(Region(subregions=[SubRegion(FeatureLocation(10, 40), "test")]))
        with self.assertRaises(ValueError):
            record.add_region(Region(subregions=[SubRegion(FeatureLocation(0, 11), "test")]))
        # ok, since ends aren't inclusive
        record.add_region(Region(subregions=[SubRegion(FeatureLocation(0, 10), "test")]))

    def test_cds_cluster_linkage(self):
        record = Record("A"*200)
        for start, end in [(50, 100), (10, 90), (0, 9), (150, 200)]:
            record.add_cds_feature(DummyCDS(start, end))
        for start, end in [(10, 120), (5, 110), (10, 160), (45, 200)]:
            record.clear_clusters()
            cluster = helpers.DummyCluster(start, end)
            record.add_cluster(cluster)
            assert len(cluster.cds_children) == 2
            for cds in cluster.cds_children:
                assert cds.overlaps_with(cluster)

    def test_orphaned_cluster_number(self):
        record = Record("A"*1000)
        cluster = helpers.DummyCluster(0, 1000)
        with self.assertRaisesRegex(ValueError, "Cluster not contained in record"):
            print(record.get_cluster_number(cluster))

    def test_orphaned_supercluster_number(self):
        record = Record("A"*1000)
        cluster = helpers.DummyCluster(0, 1000)
        supercluster = SuperCluster(SuperCluster.kinds.SINGLE, [cluster])
        with self.assertRaisesRegex(ValueError, "SuperCluster not contained in record"):
            print(record.get_supercluster_number(supercluster))

    def test_orphaned_subregion_number(self):
        record = Record(Seq("A" * 1000))
        subregion = SubRegion(FeatureLocation(0, 1000), "test")
        with self.assertRaisesRegex(ValueError, "SubRegion not contained in record"):
            print(record.get_subregion_number(subregion))

    def test_orphaned_region_number(self):
        record = Record(Seq("A" * 1000))
        subregion = SubRegion(FeatureLocation(0, 1000), "test")
        region = Region(subregions=[subregion])
        with self.assertRaisesRegex(ValueError, "Region not contained in record"):
            print(record.get_region_number(region))

# since we're about to test assigning to non-slots, shut pylint up
# pylint: disable=assigning-non-slot
    def test_membership(self):
        location = FeatureLocation(0, 3, strand=1)
        # Features don't have locus tags
        with self.assertRaises(AttributeError):
            Feature(location, feature_type="none").locus_tag = "something"
        # CDSFeatures don't have an 'other_value'
        with self.assertRaises(AttributeError):
            CDSFeature(location, translation="A", gene="a").other_value = 1
# pylint: enable=assigning-non-slot

    def test_gc_content(self):
        # pure
        for char in "ATatN":
            assert Record(Seq(char * 100)).get_gc_content() == 0.
        for char in "CGcg":
            assert Record(Seq(char * 100)).get_gc_content() == 1.

        # mixed
        self.assertAlmostEqual(Record(Seq(("A"*50) + ("G"*50))).get_gc_content(), 0.5)
        self.assertAlmostEqual(Record(Seq(("T"*25) + ("C"*75))).get_gc_content(), 0.75)

        with self.assertRaisesRegex(ValueError, "empty sequence"):
            Record().get_gc_content()
        with self.assertRaisesRegex(ValueError, "empty sequence"):
            Record("").get_gc_content()

    def test_read_from_file(self):
        # very basic testing to ensure that the file IO itself functions
        recs = Record.from_genbank(helpers.get_path_to_nisin_genbank())
        assert len(recs) == 1
        rec = recs[0]
        assert rec.get_feature_count() == 24
        assert len(rec.get_cds_features()) == 11
        assert isinstance(rec.get_cds_by_name("nisB"), CDSFeature)


class TestClusterManipulation(unittest.TestCase):
    def setUp(self):
        self.record = Record(Seq("A" * 1000))
        self.cluster = Cluster(FeatureLocation(8, 71, strand=1),
                               FeatureLocation(3, 76, strand=1), tool="test",
                               cutoff=17, neighbourhood_range=5, product='a', detection_rule="a")

    def add_cds_features(self):
        outside = DummyCDS(100, 120, locus_tag="outside")
        inside = DummyCDS(20, 40, locus_tag="inside")
        partial = DummyCDS(50, 140, locus_tag="partial")
        self.record.add_cds_feature(outside)
        self.record.add_cds_feature(inside)
        self.record.add_cds_feature(partial)
        return inside

    def test_add_cluster(self):
        assert not self.record.get_clusters()
        self.record.add_cluster(self.cluster)
        assert self.record.get_clusters() == (self.cluster,)
        assert self.record.get_cluster_number(self.cluster) == 1

    def test_add_invalid(self):
        with self.assertRaises(AssertionError):
            self.record.add_cluster(self.record)
        with self.assertRaises(AssertionError):
            self.record.add_cluster(None)
        with self.assertRaises(AssertionError):
            self.record.add_cluster(Feature(self.cluster.location, feature_type="cluster"))

    def test_add_feature(self):
        assert not self.record.get_clusters()
        self.record.add_feature(self.cluster)
        assert self.record.get_clusters() == (self.cluster,)
        assert self.record.get_cluster_number(self.cluster) == 1

    def test_clear_clusters(self):
        self.record.add_cluster(self.cluster)
        assert self.record.get_clusters()
        self.record.clear_clusters()
        assert not self.record.get_clusters()

    def test_cds_linking_cds_first(self):
        inside = self.add_cds_features()

        assert not self.cluster.cds_children
        assert self.cluster.parent_record is None
        self.record.add_cluster(self.cluster)
        assert self.cluster.parent_record is self.record
        assert self.cluster.cds_children == (inside,)

    def test_cds_linking_cluster_first(self):
        self.record.add_cluster(self.cluster)
        assert not self.cluster.cds_children

        inside = self.add_cds_features()
        assert self.cluster.cds_children == (inside,)


class TestSuperClusterManipulation(unittest.TestCase):
    def setUp(self):
        self.record = Record(Seq("A" * 1000))
        self.cluster = Cluster(FeatureLocation(8, 71), FeatureLocation(3, 76), tool="test",
                               cutoff=17, neighbourhood_range=5, product='a', detection_rule="a")
        self.record.add_cluster(self.cluster)
        self.supercluster = SuperCluster(SuperCluster.kinds.SINGLE, [self.cluster])

    def add_cds_features(self):
        outside = DummyCDS(100, 120, locus_tag="outside")
        inside = DummyCDS(20, 40, locus_tag="inside")
        partial = DummyCDS(50, 140, locus_tag="partial")
        self.record.add_cds_feature(outside)
        self.record.add_cds_feature(inside)
        self.record.add_cds_feature(partial)
        return inside

    def test_add_supercluster(self):
        assert not self.record.get_superclusters()
        self.record.add_supercluster(self.supercluster)
        assert self.record.get_superclusters() == (self.supercluster,)
        assert self.record.get_supercluster_number(self.supercluster) == 1

    def test_add_invalid(self):
        with self.assertRaises(AssertionError):
            self.record.add_supercluster(self.record)
        with self.assertRaises(AssertionError):
            self.record.add_supercluster(None)
        with self.assertRaises(AssertionError):
            self.record.add_supercluster(self.cluster)

    def test_add_feature(self):
        assert not self.record.get_superclusters()
        self.record.add_feature(self.supercluster)
        assert self.record.get_superclusters() == (self.supercluster,)
        assert self.record.get_supercluster_number(self.supercluster) == 1

    def test_clear(self):
        assert self.cluster.parent == self.supercluster
        self.record.add_supercluster(self.supercluster)
        assert self.record.get_superclusters()
        self.record.clear_superclusters()
        assert not self.record.get_superclusters()
        assert self.cluster.parent is None

    def test_creation_empty(self):
        empty_record = Record(Seq("A" * 100))
        assert not empty_record.get_superclusters()
        assert empty_record.create_superclusters() == 0
        assert not empty_record.get_superclusters()

    def test_creation_single(self):
        assert self.record.create_superclusters() == 1
        supercluster = self.cluster.parent
        assert supercluster.location == self.supercluster.location
        assert supercluster.kind == SuperCluster.kinds.SINGLE

    def test_add_biopython(self):
        bio = self.supercluster.to_biopython()[0]
        with self.assertRaisesRegex(ValueError, "cannot be directly added from biopython"):
            self.record.add_biopython_feature(bio)

    def test_cds_linking_cds_first(self):
        inside = self.add_cds_features()

        assert not self.supercluster.cds_children
        self.record.add_supercluster(self.supercluster)
        assert self.supercluster.cds_children == (inside,)

    def test_cds_linking_supercluster_first(self):
        self.record.add_supercluster(self.supercluster)
        assert not self.supercluster.cds_children

        inside = self.add_cds_features()
        assert self.supercluster.cds_children == (inside,)


class TestSubRegionManipulation(unittest.TestCase):
    def setUp(self):
        self.record = Record(Seq("A" * 1000))
        self.subregion = SubRegion(FeatureLocation(100, 200), tool="test")

    def add_cds_features(self):
        outside = DummyCDS(20, 40, locus_tag="outside")
        inside = DummyCDS(120, 140, locus_tag="inside")
        partial = DummyCDS(120, 240, locus_tag="partial")
        self.record.add_cds_feature(outside)
        self.record.add_cds_feature(inside)
        self.record.add_cds_feature(partial)
        return inside

    def test_add_subregion(self):
        assert not self.record.get_subregions()
        self.record.add_subregion(self.subregion)
        assert self.record.get_subregions() == (self.subregion,)
        assert self.record.get_subregion_number(self.subregion) == 1

    def test_add_feature(self):
        assert not self.record.get_subregions()
        self.record.add_feature(self.subregion)
        assert self.record.get_subregions() == (self.subregion,)
        assert self.record.get_subregion_number(self.subregion) == 1

    def test_add_invalid(self):
        with self.assertRaises(AssertionError):
            self.record.add_subregion(self.record)
        with self.assertRaises(AssertionError):
            self.record.add_supercluster(None)
        with self.assertRaises(AssertionError):
            self.record.add_subregion(Feature(FeatureLocation(100, 200), feature_type="subregion"))

    def test_clear(self):
        self.record.add_subregion(self.subregion)
        assert self.record.get_subregions() == (self.subregion,)
        self.record.clear_subregions()
        assert not self.record.get_subregions()

    def test_cds_linking_cds_first(self):
        inside = self.add_cds_features()

        assert not self.subregion.cds_children
        self.record.add_subregion(self.subregion)
        assert self.subregion.cds_children == (inside,)

    def test_cds_linking_subregion_first(self):
        self.record.add_subregion(self.subregion)
        assert not self.subregion.cds_children

        inside = self.add_cds_features()
        assert self.subregion.cds_children == (inside,)


class TestRegionManipulation(unittest.TestCase):
    def setUp(self):
        self.record = Record(Seq("A" * 1000))
        self.cds = DummyCDS(8, 71, locus_tag="test")
        self.record.add_cds_feature(self.cds)
        self.cluster = Cluster(FeatureLocation(8, 71), FeatureLocation(3, 76), tool="test",
                               cutoff=17, neighbourhood_range=5, product='a', detection_rule="a")
        self.record.add_cluster(self.cluster)
        self.supercluster = SuperCluster(SuperCluster.kinds.SINGLE, [self.cluster])
        self.record.add_supercluster(self.supercluster)
        self.subregion = SubRegion(FeatureLocation(200, 300), tool="test")
        self.record.add_subregion(self.subregion)
        self.region_sup = Region(superclusters=[self.supercluster])
        self.region_sub = Region(subregions=[self.subregion])
        self.region_both = Region(superclusters=[self.supercluster],
                                  subregions=[self.subregion])

    def test_add_region(self):
        assert not self.record.get_regions()
        self.record.add_region(self.region_sup)
        assert self.record.get_regions() == (self.region_sup,)
        assert self.record.get_region_number(self.region_sup) == 1

    def test_cds_linking(self):
        assert self.cds.region is None

        assert not self.cds.is_contained_by(self.region_sub)
        self.record.add_region(self.region_sub)
        assert self.cds.region is None
        new_cds = DummyCDS(220, 240, locus_tag="add_test")
        assert new_cds.region is None
        self.record.add_cds_feature(new_cds)
        assert new_cds.region is self.region_sub
        assert new_cds in self.region_sub.cds_children

        assert self.cds.is_contained_by(self.region_sup)
        self.record.add_region(self.region_sup)
        assert self.cds.region is self.region_sup
        assert self.cds in self.region_sup.cds_children

        self.record.clear_regions()
        assert self.cds.region is None

    def test_add_overlapping(self):
        assert not self.record.get_regions()
        self.record.add_region(self.region_sup)
        self.record.add_region(self.region_sub)
        with self.assertRaisesRegex(ValueError, "regions cannot overlap"):
            self.record.add_region(self.region_both)

    def test_add_feature(self):
        assert not self.record.get_regions()
        self.record.add_feature(self.region_sup)
        assert self.record.get_regions() == (self.region_sup,)
        assert self.record.get_region_number(self.region_sup) == 1

    def test_add_invalid(self):
        with self.assertRaises(AssertionError):
            self.record.add_subregion(self.record)
        with self.assertRaises(AssertionError):
            self.record.add_supercluster(None)
        with self.assertRaises(AssertionError):
            self.record.add_supercluster(self.subregion)

    def test_clear(self):
        self.record.add_region(self.region_sup)
        assert self.record.get_regions()
        self.record.clear_regions()
        assert not self.record.get_regions()

    def test_superclusters_cleared(self):
        self.record.add_region(self.region_sup)
        self.record.add_region(self.region_sub)
        assert self.record.get_superclusters()
        assert self.record.get_subregions()
        assert len(self.record.get_regions()) == 2
        self.record.clear_superclusters()
        assert len(self.record.get_regions()) == 1
        assert self.record.get_regions()[0].location == self.region_sub.location

    def test_subregions_cleared(self):
        self.record.add_region(self.region_sup)
        self.record.add_region(self.region_sub)
        assert self.record.get_superclusters()
        assert self.record.get_subregions()
        assert len(self.record.get_regions()) == 2
        self.record.clear_subregions()
        assert len(self.record.get_regions()) == 1
        assert self.record.get_regions()[0].location == self.region_sup.location

    def test_creation_empty(self):
        empty_record = Record(Seq("A" * 100))
        assert not empty_record.get_regions()
        assert empty_record.create_regions() == 0
        assert not empty_record.get_regions()

    def test_creation_independent(self):
        assert not self.record.get_regions()
        assert self.record.create_regions() == 2

        regions = self.record.get_regions()
        assert len(regions) == 2

        assert regions[0].location == self.region_sup.location
        assert regions[0].superclusters == (self.supercluster,)
        assert not regions[0].subregions

        assert regions[1].location == self.region_sub.location
        assert not regions[1].superclusters
        assert regions[1].subregions == (self.subregion,)

    def test_creation_overlapping(self):
        extra_sup = SuperCluster(SuperCluster.kinds.SINGLE, [self.cluster])
        self.record.add_supercluster(extra_sup)
        extra_sub = SubRegion(FeatureLocation(50, 250), tool="test")
        self.record.add_subregion(extra_sub)
        assert not self.record.get_regions()
        self.record.create_regions()
        assert len(self.record.get_regions()) == 1
        region = self.record.get_regions()[0]
        assert region.location == FeatureLocation(3, 300)
        assert region.superclusters == (extra_sup, self.supercluster)
        assert region.subregions == (extra_sub, self.subregion)

    def test_add_biopython(self):
        bio = self.region_sup.to_biopython()[0]
        with self.assertRaisesRegex(ValueError, "cannot be directly added from biopython"):
            self.record.add_biopython_feature(bio)
