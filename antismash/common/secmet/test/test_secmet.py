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

from antismash.common.test.helpers import get_path_to_nisin_genbank

from ..features import (
    CDSFeature,
    CandidateCluster,
    Feature,
    Region,
    Protocluster,
    SubRegion,
)
from ..errors import SecmetInvalidInputError
from .helpers import (
    DummyAntismashDomain,
    DummyCandidateCluster,
    DummyCDS,
    DummyCDSMotif,
    DummyPFAMDomain,
    DummyProtocluster,
    DummyRegion,
    DummySubRegion,
)
from ..qualifiers import SecMetQualifier, GeneFunction
from ..record import Record


class TestConversion(unittest.TestCase):
    def test_record_conversion_from_biopython(self):
        before = list(Bio.SeqIO.parse(get_path_to_nisin_genbank(), "genbank"))[0]
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
        assert type_counts["cluster"] == len(record.get_protoclusters())
        assert type_counts["aSDomain"] == len(record.get_antismash_domains())

    def test_protein_sequences_caught(self):
        before = list(Bio.SeqIO.parse(get_path_to_nisin_genbank(), "genbank"))[0]

        # as a sanity check, make sure it's a seq and it functions as expected
        assert isinstance(before.seq, Seq)
        Record.from_biopython(before, taxon="bacteria")

        before.seq = Seq("AAAA", IUPACProtein())
        with self.assertRaisesRegex(ValueError, "protein records are not supported"):
            Record.from_biopython(before, taxon="bacteria")

    def test_missing_locations_caught(self):
        rec = list(Bio.SeqIO.parse(get_path_to_nisin_genbank(), "genbank"))[0]
        Record.from_biopython(rec, taxon="bacteria")
        rec.features.append(SeqFeature(None, type="broken"))
        with self.assertRaisesRegex(SecmetInvalidInputError, "missing or invalid location"):
            Record.from_biopython(rec, taxon="bacteria")


class TestStripping(unittest.TestCase):
    def setUp(self):
        self.rec = Record(Seq("A"*20))
        self.cds = DummyCDS(locus_tag="test")
        self.rec.add_cds_feature(self.cds)

    def test_protoclusters(self):
        assert not self.rec.get_protoclusters()
        self.rec.add_protocluster(DummyProtocluster())
        assert self.rec.get_protoclusters()
        self.rec.strip_antismash_annotations()
        assert not self.rec.get_protoclusters()

    def test_candidate_clusters(self):
        assert not self.rec.get_candidate_clusters()
        self.rec.add_candidate_cluster(DummyCandidateCluster())
        assert self.rec.get_candidate_clusters()
        self.rec.strip_antismash_annotations()
        assert not self.rec.get_candidate_clusters()

    def test_subregions(self):
        assert not self.rec.get_subregions()
        self.rec.add_subregion(DummySubRegion())
        assert self.rec.get_subregions()
        self.rec.strip_antismash_annotations()
        assert not self.rec.get_subregions()

    def test_regions(self):
        assert not self.rec.get_regions()
        self.rec.add_region(DummyRegion())
        assert self.rec.get_regions()
        self.rec.strip_antismash_annotations()
        assert not self.rec.get_regions()

    def test_antismash_domains(self):
        assert not self.rec.get_antismash_domains()
        self.rec.add_antismash_domain(DummyAntismashDomain())
        assert self.rec.get_antismash_domains()
        self.rec.strip_antismash_annotations()
        assert not self.rec.get_antismash_domains()

    def test_pfams(self):
        assert not self.rec.get_pfam_domains()
        self.rec.add_pfam_domain(DummyPFAMDomain())
        assert self.rec.get_pfam_domains()
        self.rec.strip_antismash_annotations()
        assert not self.rec.get_pfam_domains()

    def test_cds_motifs(self):
        assert not self.rec.get_cds_motifs()
        motif = DummyCDSMotif()
        self.rec.add_cds_motif(motif)
        motif.created_by_antismash = False
        assert self.rec.get_cds_motifs()
        self.rec.strip_antismash_annotations()
        assert self.rec.get_cds_motifs()
        motif.created_by_antismash = True
        self.rec.strip_antismash_annotations()
        assert not self.rec.get_cds_motifs()

    def test_cds_secmet(self):
        domain = SecMetQualifier.Domain("test", 1e-2, 10., 1, "tool")
        self.cds.sec_met.add_domains([domain])
        assert self.cds.sec_met.domains == [domain]
        self.rec.strip_antismash_annotations()
        assert self.cds.sec_met.domains == []

    def test_cds_nrps_pks(self):  # pylint: disable=attribute-defined-outside-init
        class HSP:  # pylint: disable=too-few-public-methods
            def __init__(self, attrs):
                self.__dict__.update(attrs)

        raw_domain = HSP({
            "hit_id": "test",
            "query_start": 1,
            "query_end": 2,
            "evalue": 1e-5,
            "bitscore": 10,
        })
        self.cds.nrps_pks.add_domain(raw_domain, "test")
        assert self.cds.nrps_pks.domains
        self.rec.strip_antismash_annotations()
        assert not self.cds.nrps_pks.domains

    def test_cds_gene_functions(self):
        assert not self.cds.gene_functions
        self.cds.gene_functions.add(GeneFunction.OTHER, "tool", "desc")
        assert self.cds.gene_functions
        self.rec.strip_antismash_annotations()
        assert not self.cds.gene_functions


class TestRecordFeatureNumbering(unittest.TestCase):
    def setUp(self):
        self.pairs = [(50, 100), (10, 40), (700, 1000), (0, 9)]
        self.locations = [FeatureLocation(start, end) for start, end in self.pairs]
        self.record = Record(Seq("A"*1000))

    def test_protocluster_numbering(self):
        features = []
        for start, end in self.pairs:
            cluster = DummyProtocluster(start, end)
            self.record.add_protocluster(cluster)
            features.append(cluster)
        features = sorted(features)
        for i, cluster in enumerate(self.record.get_protoclusters()):
            assert cluster.get_protocluster_number() == i + 1
            assert self.record.get_protocluster(i + 1) is features[i]

    def test_candidate_cluster_numbering(self):
        features = []
        for location in self.locations:
            candidate_cluster = CandidateCluster(CandidateCluster.kinds.SINGLE,
                                                 [DummyProtocluster(location.start, location.end)])
            self.record.add_candidate_cluster(candidate_cluster)
            features.append(candidate_cluster)
        features = sorted(features)
        for i, cluster in enumerate(self.record.get_candidate_clusters()):
            assert cluster.get_candidate_cluster_number() == i + 1
            assert self.record.get_candidate_cluster(i + 1) is features[i]

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
            subregion = SubRegion(location, "test")
            self.record.add_subregion(subregion)
            region = Region(subregions=[subregion])
            self.record.add_region(region)
            features.append(region)
        features = sorted(features)
        for i, cluster in enumerate(self.record.get_regions()):
            assert cluster.get_region_number() == i + 1
            assert self.record.get_region(i + 1) is features[i]
            assert cluster.to_biopython()[0].qualifiers["region_number"] == [str(i + 1)]


class TestRecord(unittest.TestCase):
    def test_overlapping_regions(self):
        record = Record("A"*40)
        record.add_region(Region(subregions=[SubRegion(FeatureLocation(10, 40), "test")]))
        with self.assertRaises(ValueError):
            record.add_region(Region(subregions=[SubRegion(FeatureLocation(0, 11), "test")]))
        # ok, since ends aren't inclusive
        record.add_region(Region(subregions=[SubRegion(FeatureLocation(0, 10), "test")]))

    def test_cds_protocluster_linkage(self):
        record = Record("A"*200)
        for start, end in [(50, 100), (10, 90), (0, 9), (150, 200)]:
            record.add_cds_feature(DummyCDS(start, end))
        for start, end in [(10, 120), (5, 110), (10, 160), (45, 200)]:
            record.clear_protoclusters()
            cluster = DummyProtocluster(start, end)
            record.add_protocluster(cluster)
            assert len(cluster.cds_children) == 2
            for cds in cluster.cds_children:
                assert cds.overlaps_with(cluster)

    def test_orphaned_protocluster_number(self):
        record = Record("A"*1000)
        cluster = DummyProtocluster(0, 1000)
        with self.assertRaisesRegex(ValueError, "Protocluster not contained in record"):
            print(record.get_protocluster_number(cluster))

    def test_orphaned_candidate_cluster_number(self):
        record = Record("A"*1000)
        cluster = DummyProtocluster(0, 1000)
        candidate_cluster = CandidateCluster(CandidateCluster.kinds.SINGLE, [cluster])
        with self.assertRaisesRegex(ValueError, "CandidateCluster not contained in record"):
            print(record.get_candidate_cluster_number(candidate_cluster))

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
        recs = Record.from_genbank(get_path_to_nisin_genbank())
        assert len(recs) == 1
        rec = recs[0]
        assert rec.get_feature_count() == 24
        assert len(rec.get_cds_features()) == 11
        assert isinstance(rec.get_cds_by_name("nisB"), CDSFeature)

    def test_multiple_colocated_non_as_motifs(self):
        rec = Record(Seq("A" * 100))
        assert not rec.get_cds_motifs()
        domain = SeqFeature(FeatureLocation(0, 6, -1), type="CDS_motif")
        rec.add_biopython_feature(domain)
        motifs = rec.get_cds_motifs()
        assert len(motifs) == 1
        assert motifs[0].domain_id == "non_aS_motif_0_6_1"

        rec.add_biopython_feature(domain)
        motifs = rec.get_cds_motifs()
        assert len(motifs) == 2
        for i, motif in enumerate(motifs):
            assert motif.domain_id == "non_aS_motif_0_6_%s" % (i + 1)

    def test_seq_types(self):
        first = Record("A" * 20)
        assert isinstance(first.seq, Seq)
        second = Record(Seq("A" * 20))
        assert isinstance(second.seq, Seq)
        assert first.seq == second.seq


class TestCDSFetchByLocation(unittest.TestCase):
    def setUp(self):
        self.record = Record(Seq("A"*140))
        self.record.add_cds_feature(DummyCDS(10, 40, strand=1))
        self.record.add_cds_feature(DummyCDS(110, 140, strand=-1))
        self.func = self.record.get_cds_features_within_location

    def test_no_hits(self):
        assert not self.func(FeatureLocation(50, 90))

    def test_multiple_hits(self):
        assert self.func(FeatureLocation(10, 140)) == list(self.record.get_cds_features())

    def test_self_hits_same_strand(self):
        for cds in self.record.get_cds_features():
            assert self.func(cds.location, with_overlapping=False) == [cds]
            assert self.func(cds.location, with_overlapping=True) == [cds]

    def test_self_hits_opposite_strand(self):
        for cds in self.record.get_cds_features():
            loc = FeatureLocation(cds.location.start, cds.location.end, cds.location.strand * -1)
            assert self.func(loc, with_overlapping=False) == [cds]
            assert self.func(loc, with_overlapping=True) == [cds]

    def test_self_hits_no_strand(self):
        for cds in self.record.get_cds_features():
            loc = FeatureLocation(cds.location.start, cds.location.end)
            assert self.func(loc, with_overlapping=False) == [cds]
            assert self.func(loc, with_overlapping=True) == [cds]

    def test_shared_starts(self):
        starter = self.record.get_cds_features()[0]
        longer = DummyCDS(10, 80)
        shorter = DummyCDS(10, 20)
        self.record.add_cds_feature(longer)
        self.record.add_cds_feature(shorter)
        assert self.func(starter.location, with_overlapping=False) == [shorter, starter]
        assert self.func(starter.location, with_overlapping=True) == list(self.record.get_cds_features()[:3])

    def test_inception_cds(self):
        inner = self.record.get_cds_features()[1]
        outer = DummyCDS(90, 200, strand=1)
        self.record.add_cds_feature(outer)
        loc = FeatureLocation(110, 140)
        assert self.func(loc, with_overlapping=False) == [inner]


class TestClusterManipulation(unittest.TestCase):
    def setUp(self):
        self.record = Record(Seq("A" * 1000))
        self.protocluster = Protocluster(FeatureLocation(8, 71, strand=1),
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

    def test_add_protocluster(self):
        assert not self.record.get_protoclusters()
        self.record.add_protocluster(self.protocluster)
        assert self.record.get_protoclusters() == (self.protocluster,)
        assert self.record.get_protocluster_number(self.protocluster) == 1

    def test_add_invalid(self):
        with self.assertRaises(AssertionError):
            self.record.add_protocluster(self.record)
        with self.assertRaises(AssertionError):
            self.record.add_protocluster(None)
        with self.assertRaises(AssertionError):
            self.record.add_protocluster(Feature(self.protocluster.location, feature_type="protocluster"))

    def test_add_feature(self):
        assert not self.record.get_protoclusters()
        self.record.add_feature(self.protocluster)
        assert self.record.get_protoclusters() == (self.protocluster,)
        assert self.record.get_protocluster_number(self.protocluster) == 1

    def test_clear_protoclusters(self):
        self.record.add_protocluster(self.protocluster)
        assert self.record.get_protoclusters()
        self.record.clear_protoclusters()
        assert not self.record.get_protoclusters()

    def test_cds_linking_cds_first(self):
        inside = self.add_cds_features()

        assert not self.protocluster.cds_children
        assert self.protocluster.parent_record is None
        self.record.add_protocluster(self.protocluster)
        assert self.protocluster.parent_record is self.record
        assert self.protocluster.cds_children == (inside,)

    def test_cds_linking_cluster_first(self):
        self.record.add_protocluster(self.protocluster)
        assert not self.protocluster.cds_children

        inside = self.add_cds_features()
        assert self.protocluster.cds_children == (inside,)


class TestCandidateClusterManipulation(unittest.TestCase):
    def setUp(self):
        self.record = Record(Seq("A" * 1000))
        self.protocluster = Protocluster(FeatureLocation(8, 71), FeatureLocation(3, 76), tool="test",
                                         cutoff=17, neighbourhood_range=5, product='a', detection_rule="a")
        self.record.add_protocluster(self.protocluster)
        self.candidate_cluster = CandidateCluster(CandidateCluster.kinds.SINGLE, [self.protocluster])

    def add_cds_features(self):
        outside = DummyCDS(100, 120, locus_tag="outside")
        inside = DummyCDS(20, 40, locus_tag="inside")
        partial = DummyCDS(50, 140, locus_tag="partial")
        self.record.add_cds_feature(outside)
        self.record.add_cds_feature(inside)
        self.record.add_cds_feature(partial)
        return inside

    def test_add_candidate_cluster(self):
        assert not self.record.get_candidate_clusters()
        self.record.add_candidate_cluster(self.candidate_cluster)
        assert self.record.get_candidate_clusters() == (self.candidate_cluster,)
        assert self.record.get_candidate_cluster_number(self.candidate_cluster) == 1

    def test_add_invalid(self):
        with self.assertRaises(AssertionError):
            self.record.add_candidate_cluster(self.record)
        with self.assertRaises(AssertionError):
            self.record.add_candidate_cluster(None)
        with self.assertRaises(AssertionError):
            self.record.add_candidate_cluster(self.protocluster)

    def test_add_feature(self):
        assert not self.record.get_candidate_clusters()
        self.record.add_feature(self.candidate_cluster)
        assert self.record.get_candidate_clusters() == (self.candidate_cluster,)
        assert self.record.get_candidate_cluster_number(self.candidate_cluster) == 1

    def test_clear(self):
        assert self.protocluster.parent == self.candidate_cluster
        self.record.add_candidate_cluster(self.candidate_cluster)
        assert self.record.get_candidate_clusters()
        self.record.clear_candidate_clusters()
        assert not self.record.get_candidate_clusters()
        assert self.protocluster.parent is None

    def test_creation_empty(self):
        empty_record = Record(Seq("A" * 100))
        assert not empty_record.get_candidate_clusters()
        assert empty_record.create_candidate_clusters() == 0
        assert not empty_record.get_candidate_clusters()

    def test_creation_single(self):
        assert self.record.create_candidate_clusters() == 1
        candidate_cluster = self.protocluster.parent
        assert candidate_cluster.location == self.candidate_cluster.location
        assert candidate_cluster.kind == CandidateCluster.kinds.SINGLE

    def test_add_biopython(self):
        bio = self.candidate_cluster.to_biopython()[0]
        with self.assertRaisesRegex(ValueError, "cannot be directly added from biopython"):
            self.record.add_biopython_feature(bio)

    def test_cds_linking_cds_first(self):
        inside = self.add_cds_features()

        assert not self.candidate_cluster.cds_children
        self.record.add_candidate_cluster(self.candidate_cluster)
        assert self.candidate_cluster.cds_children == (inside,)

    def test_cds_linking_candidate_cluster_first(self):
        self.record.add_candidate_cluster(self.candidate_cluster)
        assert not self.candidate_cluster.cds_children

        inside = self.add_cds_features()
        assert self.candidate_cluster.cds_children == (inside,)


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
            self.record.add_candidate_cluster(None)
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
        self.protocluster = Protocluster(FeatureLocation(8, 71), FeatureLocation(3, 76), tool="test",
                                         cutoff=17, neighbourhood_range=5, product='a', detection_rule="a")
        self.record.add_protocluster(self.protocluster)
        self.candidate_cluster = CandidateCluster(CandidateCluster.kinds.SINGLE, [self.protocluster])
        self.record.add_candidate_cluster(self.candidate_cluster)
        self.subregion = SubRegion(FeatureLocation(200, 300), tool="test")
        self.record.add_subregion(self.subregion)
        self.region_sup = Region(candidate_clusters=[self.candidate_cluster])
        self.region_sub = Region(subregions=[self.subregion])
        self.region_both = Region(candidate_clusters=[self.candidate_cluster],
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
            self.record.add_candidate_cluster(None)
        with self.assertRaises(AssertionError):
            self.record.add_candidate_cluster(self.subregion)

    def test_clear(self):
        self.record.add_region(self.region_sup)
        assert self.record.get_regions()
        self.record.clear_regions()
        assert not self.record.get_regions()

    def test_candidate_clusters_cleared(self):
        self.record.add_region(self.region_sup)
        self.record.add_region(self.region_sub)
        assert self.record.get_candidate_clusters()
        assert self.record.get_subregions()
        assert len(self.record.get_regions()) == 2
        self.record.clear_candidate_clusters()
        assert len(self.record.get_regions()) == 1
        assert self.record.get_regions()[0].location == self.region_sub.location

    def test_subregions_cleared(self):
        self.record.add_region(self.region_sup)
        self.record.add_region(self.region_sub)
        assert self.record.get_candidate_clusters()
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
        assert regions[0].candidate_clusters == (self.candidate_cluster,)
        assert not regions[0].subregions

        assert regions[1].location == self.region_sub.location
        assert not regions[1].candidate_clusters
        assert regions[1].subregions == (self.subregion,)

    def test_creation_overlapping(self):
        extra_sup = CandidateCluster(CandidateCluster.kinds.SINGLE, [self.protocluster])
        self.record.add_candidate_cluster(extra_sup)
        extra_sub = SubRegion(FeatureLocation(50, 250), tool="test")
        self.record.add_subregion(extra_sub)
        assert not self.record.get_regions()
        self.record.create_regions()
        assert len(self.record.get_regions()) == 1
        region = self.record.get_regions()[0]
        assert region.location == FeatureLocation(3, 300)
        assert region.candidate_clusters == (extra_sup, self.candidate_cluster)
        assert region.subregions == (extra_sub, self.subregion)

    def test_creation_ordering(self):
        extra_cluster = Protocluster(FeatureLocation(820, 850), FeatureLocation(800, 870), tool="test",
                                     cutoff=17, neighbourhood_range=5, product='a', detection_rule="a")
        extra_sup = CandidateCluster(CandidateCluster.kinds.SINGLE, [extra_cluster])
        self.subregion.location = FeatureLocation(50, 100)
        # make sure the subregion overlaps with the first candidate cluster, but not the second
        # otherwise the test is no good
        assert self.record.get_subregion(0).overlaps_with(self.candidate_cluster)
        assert not self.record.get_subregion(0).overlaps_with(extra_sup)
        self.record.add_candidate_cluster(extra_sup)

        assert not self.record.get_regions()
        self.record.create_regions()
        assert len(self.record.get_regions()) == 2

        region = self.record.get_regions()[0]
        assert region.location == FeatureLocation(3, 100)
        assert region.candidate_clusters == (self.candidate_cluster,)
        assert region.subregions == (self.subregion,)

        region = self.record.get_regions()[1]
        assert region.location == FeatureLocation(800, 870)
        assert region.candidate_clusters == (extra_sup,)
        assert region.subregions == tuple()

    def test_add_biopython(self):
        bio = self.region_sup.to_biopython()[0]
        # it can be converted with a record, but it won't have a region number
        assert "region_number" not in bio.qualifiers
        with self.assertRaisesRegex(ValueError, "cannot be directly added from biopython"):
            self.record.add_biopython_feature(bio)


class TestCDSUniqueness(unittest.TestCase):
    def test_same_location(self):
        record = Record("A" * 100)
        cds = CDSFeature(FeatureLocation(0, 6, 1), locus_tag="test", translation="MA")
        record.add_cds_feature(cds)
        with self.assertRaisesRegex(ValueError, "same name for mapping"):
            record.add_cds_feature(cds)

        cds = CDSFeature(FeatureLocation(3, 9, 1), locus_tag="test", protein_id="prot", translation="MA")
        record.add_cds_feature(cds)
        assert cds.locus_tag == "test_prot"

        # still a duplicate
        cds = CDSFeature(FeatureLocation(3, 9, 1), locus_tag="test", protein_id="prot", translation="MA")
        with self.assertRaisesRegex(ValueError, "same name for mapping"):
            record.add_cds_feature(cds)

    def test_nonoverlapping_location(self):
        record = Record("A" * 100)
        cds = CDSFeature(FeatureLocation(0, 6, 1), locus_tag="test", translation="MA")
        record.add_cds_feature(cds)
        with self.assertRaisesRegex(ValueError, "same name for mapping"):
            record.add_cds_feature(cds)

        cds = CDSFeature(FeatureLocation(12, 18, 1), locus_tag="test", protein_id="prot", translation="MA")
        with self.assertRaisesRegex(ValueError, "same name for mapping"):
            record.add_cds_feature(cds)
