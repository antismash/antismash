# License: GNU Affero General Public License v3 or later
# A copy of GNU AGPL v3 should have been included in this software package in LICENSE.txt.

# for test files, silence irrelevant and noisy pylint warnings
# pylint: disable=use-implicit-booleaness-not-comparison,protected-access,missing-docstring

from tempfile import NamedTemporaryFile
import unittest

from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
from helperlibs.bio import seqio

from antismash.common.secmet import FeatureLocation, Record
from antismash.common.secmet.features import CandidateCluster, SubRegion, Region, Prepeptide
from antismash.common.secmet.features.region import helpers as region_helpers
from antismash.common.secmet.features.subregion import SideloadedSubRegion
from antismash.common.secmet.features.protocluster import SideloadedProtocluster
from antismash.common.secmet.locations import offset_location
from antismash.common.secmet.test.helpers import (
    DummyCandidateCluster,
    DummyCDS,
    DummyProtocluster,
    DummyRecord,
    DummyRegion,
    DummySubRegion,
)


def create_protocluster(start, end, product='a', category="A"):
    return DummyProtocluster(start, end, product=product, product_category=category)


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
        with self.assertRaisesRegex(AttributeError, "(can't set attribute|has no setter)"):
            self.region.subregions = (self.sub,)
        with self.assertRaisesRegex(AttributeError, "(can't set attribute|has no setter)"):
            self.region.candidate_clusters = (self.candidate,)
        with self.assertRaisesRegex(AttributeError, "(can't set attribute|has no setter)"):
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

    def test_product_categories(self):
        candidates = [DummyCandidateCluster([create_protocluster(0, 10, category="A")])]
        region = Region(candidate_clusters=candidates)
        assert region.product_categories == {"A"}

        clusters = []
        for i, cat in zip(range(2), ["B1", "A1"]):
            clusters.append(create_protocluster(i*10, (i+1)*10, product=cat.lower(), category=cat))
        candidate = DummyCandidateCluster(clusters)
        region = Region(candidate_clusters=[candidate])
        assert region.product_categories == {"A1", "B1"}

    def test_genbank(self):
        dummy_record = Record(Seq("A"*100))
        clusters = [DummyProtocluster(3, 20, core_start=6, core_end=12, product="prodA"),
                    DummyProtocluster(25, 41, core_start=30, core_end=36, product="prodB")]
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
        assert len(bio[0]) == len(region.location)  # ensure the sequence was trimmed

        # check the protocluster core qualifiers shifted to match the new record start
        bio_protoclusters = [feature for feature in bio[0].features if feature.type == "protocluster"]
        bio_cores = [feature for feature in bio[0].features if feature.type == "proto_core"]
        offset = -clusters[0].location.start
        expected_locations = [offset_location(cluster.core_location, offset) for cluster in clusters]
        for bio_p, location in zip(bio_protoclusters, expected_locations):
            assert bio_p.qualifiers["core_location"] == [str(location)]
        assert bio_cores
        for core, location in zip(bio_cores, expected_locations):
            assert core.location.start == location.start
            assert core.location.end == location.end

        rec = Record.from_biopython(bio[0], taxon="bacteria")
        assert len(rec.get_regions()) == 1
        new = rec.get_region(0)
        assert new.location.start == 3 - region.location.start
        assert new.location.end == 71 - region.location.start
        assert new.products == region.products

    def test_prepeptide_adjustment(self):
        dummy_record = Record(Seq("A"*400))
        subregion = DummySubRegion(start=100, end=300)
        dummy_record.add_subregion(subregion)
        region = Region(subregions=[subregion])
        dummy_record.add_region(region)

        leader_seq = "LEADER"
        core_seq = "CORE"
        tail_seq = "TAIL"
        translation = "".join([leader_seq, core_seq, tail_seq])

        dummy_prepeptide = Prepeptide(FeatureLocation(200, 200 + len(translation) * 3, 1),
                                      "dummy_class", core_seq, "locus", "dummy tool",
                                      leader=leader_seq, tail=tail_seq)
        dummy_record.add_feature(dummy_prepeptide)

        with NamedTemporaryFile(suffix=".gbk") as output:
            region.write_to_genbank(output.name)
            bio = list(seqio.parse(output.name))[0]
        assert len(bio.features) == 5  # region, subregion, leader, core, tail
        found = False
        for feature in bio.features:
            tail = feature.qualifiers.get("tail_location")
            leader = feature.qualifiers.get("leader_location")
            # only the combined feature will have the leader and tail,
            # the other two motifs created will be the features created for the tail and leader themselves
            if tail and leader:
                # the part locations should now be adjusted backwards 100 bases
                start = dummy_prepeptide.start - region.start
                end = dummy_prepeptide.end - region.start
                assert leader == [f"[{start}:{start + len(leader_seq)  * 3}](+)"]
                assert tail == [f"[{end - len(tail_seq) * 3}:{end}](+)"]
                found = True
        assert found, "prepeptide feature missing in conversion"

    def test_sideloaded(self):
        clusters = [create_protocluster(3, 20, "prodA"),
                    SideloadedProtocluster(FeatureLocation(25, 41), FeatureLocation(25, 41), "external", "prodB")]
        candidate = CandidateCluster(CandidateCluster.kinds.NEIGHBOURING, clusters)

        subregions = [SubRegion(FeatureLocation(35, 71), "test", 0.7),
                      SideloadedSubRegion(FeatureLocation(45, 61), "external")]

        region = Region(candidate_clusters=[candidate], subregions=subregions)
        sideloaded = region.get_sideloaded_areas()
        assert len(sideloaded) == 2
        assert sideloaded[0] is clusters[1]
        assert sideloaded[1] is subregions[1]


class TestHelpers(unittest.TestCase):
    def test_candidate_adjustment(self):
        bio_record = SeqRecord(seq=Seq("ACGT"))
        bio_region = SeqRecord(seq=Seq("ACGT"))
        candidates = [DummyCandidateCluster(candidate_number=i+1) for i in range(5)]
        contained_candidates = candidates[1:4]
        region = DummyRegion(candidate_clusters=contained_candidates)
        bio_region.features.extend(region.to_biopython())
        # force a fake parent record for full qualifier generation in biopython conversion
        record = DummyRecord()
        for candidate in candidates:
            candidate._parent_record = record
            for protocluster in candidate.protoclusters:
                protocluster._parent_record = record

        data = region_helpers.RegionData(
            start=0,
            end=10000,
            candidate_clusters=contained_candidates,
            subregions=[],
        )
        assert set(cand.get_candidate_cluster_number() for cand in data.candidate_clusters) == {2, 3, 4}
        bio_candidates = [cand.to_biopython()[0] for cand in candidates]
        as_quals = (int(cand.qualifiers["candidate_cluster_number"][0]) for cand in bio_candidates)
        assert set(as_quals) == {1, 2, 3, 4, 5}

        bio_record.features.extend(bio_candidates)
        bio_region.features.extend(bio_candidates[1:4])

        region_helpers._adjust_features(data, bio_region, bio_record)

        result_cands = [f for f in bio_region.features if f.type == DummyCandidateCluster.FEATURE_TYPE]
        assert len(result_cands) == 3
        print(bio_region.features)
        expected = {1, 2, 3}
        # the candidate features themselves must be adjusted
        assert set(int(cand.qualifiers["candidate_cluster_number"][0]) for cand in result_cands) == expected
        # as well as the region feature referencing those candidates
        bio_region_feature = [f for f in bio_region.features if f.type == Region.FEATURE_TYPE][0]
        assert set(map(int, bio_region_feature.qualifiers["candidate_cluster_numbers"])) == expected
