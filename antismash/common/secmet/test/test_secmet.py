# License: GNU Affero General Public License v3 or later
# A copy of GNU AGPL v3 should have been included in this software package in LICENSE.txt.

# for test files, silence irrelevant and noisy pylint warnings
# pylint: disable=use-implicit-booleaness-not-comparison,protected-access,missing-docstring

from collections import defaultdict
import unittest
from unittest.mock import patch

import Bio.SeqIO
from Bio.Seq import Seq
from Bio.SeqFeature import SeqFeature

from antismash.common.secmet import record as record_pkg
from antismash.common.secmet.locations import FeatureLocation, connect_locations
from antismash.common.test.helpers import get_path_to_nisin_genbank
from antismash.common.hmmscan_refinement import HMMResult

from ..features import (
    CDSFeature,
    CandidateCluster,
    Feature,
    Gene,
    Module,
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
    DummyRecord,
    DummyRegion,
    DummySubRegion,
)
from ..locations import (
    CompoundLocation,
)
from ..qualifiers import SecMetQualifier, GeneFunction
from ..record import ANTISMASH_SPECIFIC_TYPES, Record


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
        for i, feature in enumerate(before.features):
            assert id(feature) != id(after.features[i])
        for bef, aft in zip(before_features, sorted(map(str, before.features))):
            assert bef == aft

        # ensure that the counts of each match
        assert type_counts["CDS"] == len(record.get_cds_features())
        assert type_counts["PFAM_domain"] == len(record.get_pfam_domains())
        assert type_counts["cluster"] == len(record.get_protoclusters())
        assert type_counts["aSDomain"] == len(record.get_antismash_domains())

    def test_conversion_with_kwargs(self):
        bio = Record(seq=Seq("ACGT")).to_biopython()
        other_values = {"first": 5., "second": "value"}
        with patch.object(record_pkg.SeqRecord, "__init__", return_value=None) as patched_bio:
            Record.from_biopython(bio, "taxon", **other_values)
            patched_bio.assert_called_once_with(bio.seq, **other_values)

    def test_cross_origin_feature_in_linear(self):
        # valid locations
        for location in [
            FeatureLocation(0, 100, 1),
            FeatureLocation(0, 100, -1),
            CompoundLocation([FeatureLocation(0, 10, 1), FeatureLocation(90, 100, 1)]),
            CompoundLocation([FeatureLocation(90, 100, -1), FeatureLocation(0, 10, -1)]),
        ]:
            bio = Record(seq=Seq("A" * 100)).to_biopython()
            bio.features.append(SeqFeature(location, type="test"))
            Record.from_biopython(bio, taxon="bacteria")

        # invalid locations
        for location in [
            CompoundLocation([FeatureLocation(0, 10, -1), FeatureLocation(90, 100, -1)]),
            CompoundLocation([FeatureLocation(90, 100, 1), FeatureLocation(0, 10, 1)]),
        ]:
            bio = Record(seq=Seq("A" * 100)).to_biopython()
            bio.features.append(SeqFeature(location, type="test"))
            with self.assertRaisesRegex(SecmetInvalidInputError, "origin spanning .* in a linear record"):
                Record.from_biopython(bio, taxon="bacteria")

    def test_protein_sequences_caught(self):
        before = list(Bio.SeqIO.parse(get_path_to_nisin_genbank(), "genbank"))[0]

        # as a sanity check, make sure it's a seq and it functions as expected
        assert isinstance(before.seq, Seq)
        Record.from_biopython(before, taxon="bacteria")

        before.annotations["molecule_type"] = "protein"
        with self.assertRaisesRegex(ValueError, "protein records are not supported"):
            Record.from_biopython(before, taxon="bacteria")

    def test_dna_casing(self):
        before = list(Bio.SeqIO.parse(get_path_to_nisin_genbank(), "genbank"))[0]
        for molecule in ["DNA", "dna", "Dna"]:
            before.annotations["molecule_type"] = molecule
            Record.from_biopython(before, taxon="bacteria")

            before.annotations["molecule_type"] = molecule + "x"
            with self.assertRaisesRegex(ValueError, "records are not supported"):
                Record.from_biopython(before, taxon="bacteria")

    def test_missing_locations_caught(self):
        rec = list(Bio.SeqIO.parse(get_path_to_nisin_genbank(), "genbank"))[0]
        Record.from_biopython(rec, taxon="bacteria")
        rec.features.append(SeqFeature(None, type="broken"))
        with self.assertRaisesRegex(SecmetInvalidInputError, "missing or invalid location"):
            Record.from_biopython(rec, taxon="bacteria")

    def test_discard(self):
        bio = list(Bio.SeqIO.parse(get_path_to_nisin_genbank(), "genbank"))[0]
        for feature_type in ANTISMASH_SPECIFIC_TYPES:
            bio.features = [
                SeqFeature(FeatureLocation(10, 70, 1), type=feature_type),  # invalid
                SeqFeature(FeatureLocation(100, 160, 1), type=CDSFeature.FEATURE_TYPE),  # valid/non-antismash
            ]
            # without the discard flag, this needs to raise an error
            with self.assertRaisesRegex(SecmetInvalidInputError,
                                        "(missing expected qualifier: '[^'])|(requires at least one)"):
                Record.from_biopython(bio, taxon="bacteria")
            # with the discard flag, the bad antismash-specific feature should just disappear
            rec = Record.from_biopython(bio, taxon="bacteria", discard_antismash_features=True)
            assert len(list(rec.all_features)) == 1
            # and that single feature had better be the CDS
            assert len(rec.get_cds_features()) == 1


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
        self.rec.add_candidate_cluster(DummyCandidateCluster(end=len(self.rec)))
        assert self.rec.get_candidate_clusters()
        self.rec.strip_antismash_annotations()
        assert not self.rec.get_candidate_clusters()

    def test_subregions(self):
        assert not self.rec.get_subregions()
        self.rec.add_subregion(DummySubRegion(end=len(self.rec)))
        assert self.rec.get_subregions()
        self.rec.strip_antismash_annotations()
        assert not self.rec.get_subregions()

    def test_regions(self):
        assert not self.rec.get_regions()
        self.rec.add_region(DummyRegion(end=len(self.rec)))
        assert self.rec.get_regions()
        self.rec.strip_antismash_annotations()
        assert not self.rec.get_regions()

    def test_antismash_domains(self):
        assert not self.rec.get_antismash_domains()
        self.rec.add_antismash_domain(DummyAntismashDomain(tool="test"))
        assert self.rec.get_antismash_domains()
        assert self.rec.get_antismash_domains_by_tool("test")
        self.rec.strip_antismash_annotations()
        assert not self.rec.get_antismash_domains()
        assert not self.rec.get_antismash_domains_by_tool("test")

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
        raw_domain = HMMResult("test", 1, 2, 1e-5, 10)
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

    def test_antismash_domains_by_tool(self):
        assert not self.rec.get_antismash_domains()
        assert self.rec.get_antismash_domains_by_tool("a") == tuple()
        assert self.rec.get_antismash_domains_by_tool("b") == tuple()
        a = DummyAntismashDomain(tool="a")
        b = DummyAntismashDomain(tool="b")
        z = DummyAntismashDomain(tool="a")
        for i in [a, b, z]:
            self.rec.add_antismash_domain(i)
        assert self.rec.get_antismash_domains() == (a, b, z)
        assert self.rec.get_antismash_domains_by_tool("a") == (a, z)
        assert self.rec.get_antismash_domains_by_tool("b") == (b,)


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

    def test_cds_caching(self):
        record = Record("A" * 100)
        record.add_cds_feature(DummyCDS(10, 40, strand=1))
        original = record.get_cds_features()
        assert len(original) == 1
        assert record.get_cds_features() is original  # same object, since no change
        record.add_cds_feature(DummyCDS(110, 140, strand=-1))
        updated = record.get_cds_features()
        assert len(updated) == 2
        assert updated is not original

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
            assert motif.domain_id == f"non_aS_motif_0_6_{i + 1}"

    def test_seq_types(self):
        first = Record("A" * 20)
        assert isinstance(first.seq, Seq)
        second = Record(Seq("A" * 20))
        assert isinstance(second.seq, Seq)
        assert first.seq == second.seq

    def test_gc_caching(self):
        rec = Record("ATCG" * 20)
        assert rec._gc_content == -1  # cache should not be set
        # ensure 'Counter' is actually the mechanism for calculation, otherwise the test will be inaccurate
        with patch.object(record_pkg, "Counter") as patched:
            rec.get_gc_content()
            patched.assert_called_once_with(rec.seq)
        # since that will have cached the mock, reset it
        rec._gc_content = -1
        # now start the real calculation
        gc_content = rec.get_gc_content()
        self.assertAlmostEqual(gc_content, 0.5)
        assert rec._gc_content == gc_content  # cache should be updated
        # and ensure that it's cached and still the right value
        with patch.object(record_pkg, "Counter") as patched:
            assert rec.get_gc_content() == gc_content  # value should not change
            assert not patched.called  # and the calculation shouldn't have run again

    def test_gc_setting(self):
        rec = Record("A", gc_content=0.3)
        # the arg is trusted, even if it's wrong
        assert rec._gc_content == 0.3
        # and since it's already set, the getter should return it
        assert rec.get_gc_content() == 0.3

    def test_distance_between_features(self):
        record = Record("A" * 100)
        assert not record.is_circular()
        low = Feature(FeatureLocation(10, 20, 1), "test")
        high = Feature(FeatureLocation(70, 80, 1), "test")

        with patch.object(record_pkg, "get_distance_between_locations", return_value="dummy") as patched:
            # make sure what's returned is the location function's result
            assert record.get_distance_between_features(high, low) == "dummy"
            # and make sure it's called with the right args
            patched.assert_called_once_with(high.location, low.location)

        # and check with a circular record
        record._record.annotations = {"topology": "circular"}
        assert record.is_circular()

        with patch.object(record_pkg, "get_distance_between_locations", return_value="dummy") as patched:
            # again, make sure what's returned is the location function's result
            assert record.get_distance_between_features(high, low) == "dummy"
            # and make sure that the wrap point is present
            patched.assert_called_once_with(high.location, low.location, wrap_point=len(record))


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

    def test_compound_with_overlapping(self):
        location = CompoundLocation([
            FeatureLocation(20, 40, 1),
            FeatureLocation(60, 80, 1),
        ])
        for cds, overlap_expected in zip(self.record.get_cds_features(), [True, False]):
            assert not cds.is_contained_by(location)
            assert cds.overlaps_with(location) == overlap_expected

        assert self.func(location, with_overlapping=False) == []
        assert self.func(location, with_overlapping=True) == [self.record.get_cds_features()[0]]


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
        self.record.add_biopython_feature(bio)
        assert self.record.get_candidate_clusters()

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

    def test_creation_awkward(self):
        # occurred in NZ_CP017242.1
        record_length = 430631
        cross_origin = DummyCandidateCluster(clusters=[
            DummyProtocluster(core_location=CompoundLocation([
                FeatureLocation(415982, record_length, 1),
                FeatureLocation(0, 5966, 1),
            ]), neighbourhood_range=0, record_length=record_length),
        ])
        single = DummyCandidateCluster(start=8237, end=74008)
        others = [
            DummyCandidateCluster(start=188600, end=238069),
            DummyCandidateCluster(start=188600, end=231636),
            DummyCandidateCluster(start=192750, end=238069),
        ]
        record = DummyRecord(length=record_length, circular=True)
        regions = record.create_regions([cross_origin, single] + others)
        assert regions == 3

        independent, merged, cross = sorted(record.get_regions(), key=lambda x: x.start)

        assert cross.location == cross_origin.location
        assert independent.location == single.location
        assert merged.location == record.connect_locations([other.location for other in others])

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

    def test_circular_without_cross_origin(self):
        record = Record(Seq("A" * 100))
        record.annotations["topology"] = "circular"
        assert record.is_circular()

        early = SubRegion(FeatureLocation(13, 26, 1), tool="test")
        late = SubRegion(FeatureLocation(75, 97, 1), tool="test")
        assert not early.overlaps_with(late)

        for sub in [early, late]:
            record.add_subregion(sub)

        record.create_regions()
        regions = record.get_regions()
        assert len(regions) == 2
        assert regions[0].location == early.location
        assert regions[1].location == late.location

    def test_circular_creation_with_cross_origin(self):
        record = Record(Seq("A" * 100))
        record.annotations["topology"] = "circular"
        assert record.is_circular()

        early = SubRegion(FeatureLocation(13, 26, 1), tool="test")
        late = SubRegion(FeatureLocation(75, 97, 1), tool="test")
        origin = SubRegion(CompoundLocation([
            FeatureLocation(85, 100, 1),
            FeatureLocation(0, 5, 1),
        ]), tool="test")

        assert not late.overlaps_with(early)
        assert origin.overlaps_with(late)
        assert not origin.overlaps_with(early)

        for sub in [early, late, origin]:
            record.add_subregion(sub)

        record.create_regions()
        regions = record.get_regions()
        assert len(regions) == 2
        cross_origin, other = sorted(regions, key=lambda x: x.location.start)

        assert cross_origin.crosses_origin()
        assert cross_origin.start == late.start and cross_origin.end == origin.end

        assert not other.crosses_origin()
        assert other.location == early.location

    def test_creation_overlapping(self):
        extra_sup = CandidateCluster(CandidateCluster.kinds.SINGLE, [self.protocluster])
        self.record.add_candidate_cluster(extra_sup)
        extra_sub = SubRegion(FeatureLocation(50, 250), tool="test")
        self.record.add_subregion(extra_sub)
        assert not self.record.get_regions()
        self.record.create_regions()
        assert len(self.record.get_regions()) == 1
        region = self.record.get_regions()[0]
        assert region.location == FeatureLocation(3, 300, 1)
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
        assert region.location == FeatureLocation(3, 100, 1)
        assert region.candidate_clusters == (self.candidate_cluster,)
        assert region.subregions == (self.subregion,)

        region = self.record.get_regions()[1]
        assert region.location == FeatureLocation(800, 870, 1)
        assert region.candidate_clusters == (extra_sup,)
        assert region.subregions == tuple()

    def test_add_biopython(self):
        bio = self.region_sup.to_biopython()[0]
        # it can be converted with a record, but it won't have a region number
        assert "region_number" not in bio.qualifiers
        self.record.add_biopython_feature(bio)
        assert self.record.get_regions()


class TestModuleManipulation(unittest.TestCase):
    def setUp(self):
        self.record = Record(Seq("A" * 1000))
        self.cdses = [DummyCDS(8, 71), DummyCDS(80, 110), DummyCDS(100, 180)]
        self.domains = [DummyAntismashDomain(start=40, end=61),
                        DummyAntismashDomain(start=90, end=99)]
        self.domains[0].locus_tag = self.cdses[0].get_name()
        self.domains[1].locus_tag = self.cdses[1].get_name()

        for feature in self.cdses + self.domains:
            self.record.add_feature(feature)

        location = connect_locations([dom.location for dom in self.domains])
        self.module = Module(location, domains=self.domains)

    def test_add_biopython(self):
        bio = self.module.to_biopython()[0]
        self.record.add_biopython_feature(bio)
        assert self.record.get_modules()
        for cds in self.cdses[:2]:
            assert cds.modules
        assert not self.cdses[2].modules


class TestCDSUniqueness(unittest.TestCase):
    def test_same_location(self):
        record = Record("A" * 100)
        cds = CDSFeature(FeatureLocation(0, 6, 1), locus_tag="test", translation="MA")
        record.add_cds_feature(cds)
        with self.assertRaisesRegex(ValueError, "same location"):
            record.add_cds_feature(cds)

        cds = CDSFeature(FeatureLocation(3, 9, 1), locus_tag="test", protein_id="prot", translation="MA")
        record.add_cds_feature(cds)
        assert cds.locus_tag == "test_c6ced428"

        # still a duplicate
        cds = CDSFeature(FeatureLocation(3, 9, 1), locus_tag="test", protein_id="prot", translation="MA")
        with self.assertRaisesRegex(ValueError, "same location"):
            record.add_cds_feature(cds)

        # reverse strand, same coordinates are fine
        cds = CDSFeature(FeatureLocation(0, 6, -1), locus_tag="test_reverse", translation="MA")
        record.add_cds_feature(cds)

        # no locus tag, conflicting protein_id
        cds = CDSFeature(FeatureLocation(12, 18, 1), protein_id="test", translation="MA")
        with self.assertRaisesRegex(ValueError, "same name for mapping"):
            record.add_cds_feature(cds)

    def test_nonoverlapping_location(self):
        record = Record("A" * 100)
        cds = CDSFeature(FeatureLocation(0, 6, 1), locus_tag="test", translation="MA")
        record.add_cds_feature(cds)
        with self.assertRaisesRegex(ValueError, "same location"):
            record.add_cds_feature(cds)

        cds = CDSFeature(FeatureLocation(12, 18, 1), locus_tag="test", protein_id="prot", translation="MA")
        with self.assertRaisesRegex(ValueError, "same name for mapping"):
            record.add_cds_feature(cds)

        # Allowed special case: multiple CDS splice variants from the same gene
        # can exist even if they don't overlap
        gene = Gene(FeatureLocation(20, 42, 1), locus_tag="test_gene")
        record.add_gene(gene)

        splice1 = CDSFeature(FeatureLocation(20, 26, 1), locus_tag="test_gene",
                             protein_id="splice1", translation="MA")
        splice2 = CDSFeature(FeatureLocation(27, 33, 1), locus_tag="test_gene",
                             protein_id="splice2", translation="MA")
        record.add_cds_feature(splice1)
        record.add_cds_feature(splice2)

        assert splice1.locus_tag == "test_gene"
        assert splice2.locus_tag == "test_gene_7e364f92"

        with self.assertRaisesRegex(ValueError, "same location"):
            splice1 = CDSFeature(FeatureLocation(20, 26, 1), locus_tag="test_gene",
                                 protein_id="splice1", translation="MA")
            record.add_cds_feature(splice1)

        with self.assertRaisesRegex(ValueError, "same location"):
            splice2 = CDSFeature(FeatureLocation(27, 33, 1), locus_tag="test_gene",
                                 protein_id="splice2", translation="MA")
            record.add_cds_feature(splice2)

        splice3 = CDSFeature(FeatureLocation(36, 42, 1), locus_tag="test_gene", protein_id="splice2", translation="MA")
        record.add_cds_feature(splice3)
        assert splice3.locus_tag == "test_gene_85824247"


class TestIsNuclSeq(unittest.TestCase):
    def test_seq(self):
        # > 20%
        for seq in ["AGTC", "AGCTFC", "agtcfc", "AGTCFCT"]:
            assert Record.is_nucleotide_sequence(Seq(seq))
        # edge case == 20% should be failure
        assert not Record.is_nucleotide_sequence(Seq("AGFTC"))
        # and less than 20%
        assert not Record.is_nucleotide_sequence(Seq("AGFTCF"))

    def test_str(self):
        for seq in ["AGTC", "AGCTFC", "agtcfc", "AGTCFCT"]:
            assert Record.is_nucleotide_sequence(seq)
        assert not Record.is_nucleotide_sequence("AGFTC")
        assert not Record.is_nucleotide_sequence("AGFTCF")


def test_naming():
    short_name = "some_id"
    long_name = "some_really_very_long_id_needing_truncation"
    record = Record(Seq("AAA"))
    assert record.original_id is None

    assert not record.has_name("other")

    record.id = short_name
    assert record.original_id is None
    assert not record.has_name(None)
    assert not record.has_name("")

    record.original_id = long_name
    assert record.has_name(short_name)
    assert record.has_name(long_name)

    record.id = "other"
    assert record.has_name(long_name)


class TestExtension(unittest.TestCase):
    def setUp(self):
        self.record = Record("A" * 100)
        assert not self.record.is_circular()

    def set_circular(self):
        self.record._record.annotations = {"topology": "circular"}
        assert self.record.is_circular()

    def test_extending_outside_both_ends(self):
        self.set_circular()
        full_location = FeatureLocation(0, 100, 1)
        assert len(full_location) == len(self.record)

        new = self.record.extend_location(full_location, 10)
        assert new.parts is not full_location.parts  # to avoid mutability issues
        assert new.parts == full_location.parts  # but they ought to be identical

        location = FeatureLocation(2, 97, 1)
        new = self.record.extend_location(location, 10)
        assert new.parts is not location.parts
        assert new.parts == full_location.parts

    def test_simple_unbounded(self):
        for strand in [-1, 0, 1]:
            location = FeatureLocation(30, 50, strand)
            for distance in [10, 15]:
                new = self.record.extend_location(location, distance)
                assert isinstance(new, FeatureLocation)
                assert new.start == location.start - distance
                assert new.end == location.end + distance
                assert new.strand == location.strand
                # ensure the original wasn't changed
                assert new.parts is not location.parts

    def test_simple_underflow(self):
        location = FeatureLocation(10, 30, 1)
        distance = 20
        new = self.record.extend_location(location, distance)
        assert isinstance(new, FeatureLocation)
        assert new.start == 0
        assert new.end == location.end + distance
        assert new.strand == location.strand

        self.set_circular()
        assert self.record.is_circular()

        new = self.record.extend_location(location, distance)
        assert isinstance(new, CompoundLocation)
        assert new.start == 0
        assert new.end == len(self.record)
        assert new.strand == location.strand
        assert new.parts == [
            FeatureLocation((location.start - distance) % len(self.record), new.end, 1),
            FeatureLocation(0, location.end + distance, 1)
        ]
        # and ensure the original wasn't changed
        assert new.parts is not location.parts

    def test_simple_overflow(self):
        location = FeatureLocation(60, 90, 1)
        distance = 20
        new = self.record.extend_location(location, distance)
        assert isinstance(new, FeatureLocation)
        assert new.start == location.start - distance
        assert new.end == len(self.record)
        assert new.strand == location.strand

        self.set_circular()
        assert self.record.is_circular()

        new = self.record.extend_location(location, distance)
        assert isinstance(new, CompoundLocation)
        assert new.start == 0
        assert new.end == len(self.record)
        assert new.strand == location.strand
        assert new.parts == [
            FeatureLocation(location.start - distance, new.end, 1),
            FeatureLocation(0, (location.end + distance) % len(self.record), 1)
        ]

    def test_simple_too_long(self):
        location = FeatureLocation(40, 60, 1)
        distance = 200  # longer than the record itself should be fine
        new = self.record.extend_location(location, distance)
        assert new.start == 0
        assert new.end == 100

    def test_compound_unbounded(self):
        for strand in [-1, 0, 1]:
            initial_parts = [
                FeatureLocation(20, 40, strand),
                FeatureLocation(60, 80, strand),
            ]
            if strand == -1:
                location = CompoundLocation(initial_parts[::-1])
            else:
                location = CompoundLocation(initial_parts)
            for distance in [10, 15]:
                new = self.record.extend_location(location, distance, connect_result=False)
                assert isinstance(new, CompoundLocation)
                assert new.start == location.start - distance
                assert new.end == location.end + distance
                assert new.strand == location.strand
                assert len(new.parts) == len(location.parts)
                parts = sorted(location.parts, key=lambda x: x.start)
                assert parts[0].end == parts[0].end
                assert parts[1].start == parts[1].start
                # ensure the original wasn't changed
                assert parts is not location.parts

    def test_compound_underflow(self):
        location = CompoundLocation([
            FeatureLocation(10, 30, 1),
            FeatureLocation(40, 60, 1),
        ])
        distance = 20
        new = self.record.extend_location(location, distance, connect_result=False)
        assert isinstance(new, CompoundLocation)
        assert new.start == 0
        assert new.end == location.end + distance
        assert new.strand == location.strand
        assert new.parts is not location.parts
        assert new.parts[0].end == location.parts[0].end
        assert new.parts[1].start == location.parts[1].start

        self.set_circular()
        assert self.record.is_circular()

        new = self.record.extend_location(location, distance, connect_result=True)
        assert isinstance(new, CompoundLocation)
        assert new.start == 0
        assert new.end == len(self.record)
        assert new.strand == location.strand
        assert new.parts == [
            FeatureLocation(90, len(self.record), 1),
            FeatureLocation(0, location.end + distance, 1),
        ]
        assert new.parts is not location.parts

    def test_compound_overflow(self):
        location = CompoundLocation([
            FeatureLocation(40, 60, 1),
            FeatureLocation(70, 90, 1),
        ])
        distance = 20
        new = self.record.extend_location(location, distance, connect_result=False)
        assert isinstance(new, CompoundLocation)
        assert new.start == location.start - distance
        assert new.end == len(self.record)
        assert new.strand == location.strand
        assert new.parts is not location.parts
        assert new.parts[0].end == location.parts[0].end
        assert new.parts[1].start == location.parts[1].start

        self.set_circular()
        assert self.record.is_circular()

        new = self.record.extend_location(location, distance, connect_result=True)
        assert isinstance(new, CompoundLocation)
        assert new.start == 0
        assert new.end == len(self.record)
        assert new.strand == location.strand
        assert new.parts == [
            FeatureLocation(20, len(self.record), 1),
            FeatureLocation(0, 10, 1),
        ]
        assert new.parts is not location.parts

    def test_compound_too_long_linear(self):
        location = CompoundLocation([
            FeatureLocation(20, 30, 1),
            FeatureLocation(40, 50, 1),
            FeatureLocation(60, 70, 1),
        ])
        distance = 200  # longer than the record
        new = self.record.extend_location(location, distance, connect_result=False)
        # while they should extend to both edges, they shouldn't merge because they don't wrap
        assert new.parts[0].start == 0 and new.parts[0].end == 30
        assert new.parts[1] == location.parts[1]
        assert new.parts[2].start == 60 and new.parts[2].end == 100

    def test_compound_too_long_circular(self):
        location = CompoundLocation([
            FeatureLocation(20, 30, 1),
            FeatureLocation(40, 50, 1),
            FeatureLocation(60, 70, 1),
        ])
        distance = 100  # long enough for the extension to wrap and cover the disjoint areas in the middle
        self.set_circular()
        new = self.record.extend_location(location, distance)
        assert new.start == 0
        assert new.end == 100
        # and since it's connecting parts now, they should be joined
        assert len(new.parts) == 1

    def test_real(self):
        record = DummyRecord(seq="A"*212)
        location = CompoundLocation([
            FeatureLocation(112, 212, 1),
            FeatureLocation(0, 105, 1),
        ])
        result = record.extend_location(location, 100)
        assert result == FeatureLocation(0, 212, 1)


class TestExtensionConnect(unittest.TestCase):
    # these tests are just for connection, another class tests all the extensions over the origin and such
    def setUp(self):
        self.record = DummyRecord(length=100, circular=True)

    def compare(self, original, expected, distance, connect=True):
        new = self.record.extend_location(original, distance, connect_result=connect)
        assert original.strand == expected[0][-1]
        assert new.strand == original.strand
        assert len(new.parts) == len(expected)
        for i, part in enumerate(new.parts):
            start, end, strand = expected[i]
            assert part.start == start
            assert part.end == end
            assert part.strand == strand
        if len(original.parts) > 1 and len(new.parts) > 1:
            assert original.operator == new.operator

    def test_unbounded_compound(self):
        for strand in [1, -1]:
            parts = [FeatureLocation(20, 30, strand), FeatureLocation(50, 60, strand)]
            if strand == -1:
                parts.reverse()
            location = CompoundLocation(parts)
            distance = 10
            self.compare(location, [(10, 30, strand), (50, 70, strand)][::strand], distance, connect=False)
            self.compare(location, [(10, 70, strand)][::strand], distance, connect=True)

    def test_unbounded_simple(self):
        for strand in [1, -1]:
            location = FeatureLocation(20, 60, strand)
            self.compare(location, [(10, 70, strand)], distance=10, connect=False)
            self.compare(location, [(10, 70, strand)], distance=10, connect=True)

    def test_cross_origin_reverse(self):
        parts = [FeatureLocation(10, 20, -1), FeatureLocation(80, 90, -1)]
        for operator in ["join", "order"]:
            location = CompoundLocation(parts, operator=operator)
            self.compare(location, [(10, 25, -1), (75, 90, -1)], distance=5, connect=False)
            self.compare(location, [(0, 25, -1), (75, 100, -1)], distance=5, connect=True)

    def test_cross_origin_forward(self):
        parts = [FeatureLocation(80, 90, 1), FeatureLocation(10, 20, 1)]
        for operator in ["join", "order"]:
            location = CompoundLocation(parts, operator=operator)
            self.compare(location, [(75, 90, 1), (10, 25, 1)], distance=5, connect=False)
            self.compare(location, [(75, 100, 1), (0, 25, 1)], distance=5, connect=True)
