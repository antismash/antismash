# License: GNU Affero General Public License v3 or later
# A copy of GNU AGPL v3 should have been included in this software package in LICENSE.txt.

# for test files, silence irrelevant and noisy pylint warnings
# pylint: disable=use-implicit-booleaness-not-comparison,protected-access,missing-docstring

import json
import unittest

from Bio.Seq import Seq
from Bio.SeqFeature import FeatureLocation

from antismash import get_all_modules
from antismash.common.secmet.test.helpers import DummyRegion, DummySubRegion
from antismash.common.test.helpers import (
    DummyCandidateCluster,
    DummyCDS,
    DummyProtocluster,
    DummyRecord,
)
from antismash.config import build_config, destroy_config, get_config, update_config
from antismash.modules.tfbs_finder.tfbs_finder import (
    Confidence,
    Hit,
    Matrix,
    TFBSHit,
    TFBSFinderResults,
    get_valid_areas,
    load_matrices,
    filter_hits,
    get_sequence_gc_content,
    get_bg_distribution,
    PWM_PATH,
)
from antismash.modules.tfbs_finder.html_output import (
    add_neighbouring_genes,
    generate_javascript_data,
    get_sequence_matches,
)


class TestResults(unittest.TestCase):
    def setUp(self):
        self.record = DummyRecord()
        self.record.id = 'test_record'
        self.hits_by_region = {1: [TFBSHit('Test1', 1, 'TestHit', 'TGA', Confidence.STRONG, 1, 10, 10)]}

        build_config([
            "--tfbs",
            "--tfbs-pvalue", "0.0001",
            "--tfbs-range", "50",
        ], isolated=True, modules=get_all_modules())

    def tearDown(self):
        destroy_config()

    def create_results(self, record_id="test_record", pval=0.0001, start_overlap=50,
                       hits_by_region=None):
        hits_by_region = hits_by_region or self.hits_by_region

        return TFBSFinderResults(record_id, pval, start_overlap, hits_by_region)

    def test_init(self):
        results = TFBSFinderResults("test_record", 0.0001, 50, self.hits_by_region)
        assert results.record_id == self.record.id
        assert results.pvalue == 0.0001
        assert results.start_overlap == 50
        assert results.hits_by_region is self.hits_by_region

    def test_json_conversion(self):
        results = self.create_results()
        data = json.loads(json.dumps(results.to_json()))

        regenerated = TFBSFinderResults.from_json(data, self.record)
        assert regenerated.hits_by_region == results.hits_by_region
        assert regenerated.pvalue == results.pvalue
        assert regenerated.start_overlap == results.start_overlap
        assert regenerated.record_id == results.record_id

    def test_from_json_wrong_record(self):
        data = self.create_results().to_json()
        data['record_id'] = 'wrong_record'
        assert not TFBSFinderResults.from_json(data, self.record)

    def test_from_json_wrong_schema_version(self):
        data = self.create_results().to_json()
        data['schema_version'] = 'wrong_schema'
        assert not TFBSFinderResults.from_json(data, self.record)

    def test_different_pvalue(self):
        data = self.create_results().to_json()
        assert get_config().tfbs_pvalue == 0.0001
        update_config({"tfbs_pvalue": 0.01})
        assert not TFBSFinderResults.from_json(data, self.record)

    def test_different_range(self):
        data = self.create_results().to_json()
        assert get_config().tfbs_range == 50
        update_config({"tfbs_range": 60})
        assert not TFBSFinderResults.from_json(data, self.record)

    def test_add_to_incorrect_record(self):
        results = self.create_results(record_id=self.record.id)
        with self.assertRaisesRegex(ValueError,
                                    "Record to store in and record analysed don't match"):
            other = DummyRecord()
            other.id = self.record.id * 2
            results.add_to_record(other)

    def test_new_feature_from_hits(self):
        data = self.create_results(hits_by_region=self.hits_by_region)
        assert len(data.features) == 1
        feature = data.features[0]
        assert feature.strand == 1
        assert feature.location == FeatureLocation(1, 4, 1)
        assert feature.type == "misc_feature"
        assert feature.created_by_antismash


def make_dummy_matrix(name="name"):
    nuc_a = [1.0, 1.0, -1.0, -1.0, -1.0, -1.0]
    nuc_c = [-1.0, -1.0, 1.0, -1.0, -1.0, -1.0]
    nuc_g = [-1.0, -1.0, -1.0, 1.0, -1.0, -1.0]
    nuc_t = [-1.0, -1.0, -1.0, -1.0, 1.0, 1.0]
    matrix = [nuc_a, nuc_c, nuc_g, nuc_t]
    matrix_max_score = 6.0
    matrix_min_score = 3.0
    description = 'desc'
    consensus = 'AACGTT'
    return Matrix(name, matrix, matrix_min_score, matrix_max_score, description, consensus)


class TestTFBSFinder(unittest.TestCase):
    def setUp(self):
        self.hit_sequence = Seq("GAATAACGTTAGGCTCAACGTTGCT")
        self.non_hit_sequence = Seq("GTAGGGTACCCATACAGGTCCTATGCT")
        self.record = self.make_dummy_record()

    def make_dummy_record(self):
        protocluster = DummyProtocluster(start=0, end=25)
        cds1 = DummyCDS(start=0, end=13, strand=-1, locus_tag='a')
        cds2 = DummyCDS(start=23, end=25, strand=1, locus_tag='b')
        candidatecluster = DummyCandidateCluster(clusters=[protocluster])
        hit_region = DummyRegion(candidate_clusters=[candidatecluster])
        record = DummyRecord(seq=self.hit_sequence,
                             features=[cds1, cds2, protocluster, hit_region])
        return record

    def test_get_sequence_gc_content(self):
        gc_perc = get_sequence_gc_content(str(self.hit_sequence))
        assert gc_perc == 0.44

    def test_get_bg_distribution(self):
        background = get_bg_distribution(self.hit_sequence)
        assert background == (0.28, 0.22, 0.22, 0.28)
        assert get_bg_distribution(Seq("ATAA")) == (0.5, 0.0, 0.0, 0.5)
        assert get_bg_distribution(Seq("ATGC")) == (0.25, 0.25, 0.25, 0.25)
        assert get_bg_distribution(Seq("GCGC")) == (0.0, 0.5, 0.5, 0.0)

    def test_filtering(self):
        matrices = [make_dummy_matrix(name="A"), make_dummy_matrix(name="B")]
        hits = [
            [  # A
                Hit(pos=160, score=0, strand=1),
                Hit(pos=120, score=0, strand=1),
            ],
            [  # B
                Hit(pos=400, score=0, strand=-1),
                Hit(pos=300, score=0, strand=-1),
            ],
        ]
        filtered = filter_hits(matrices, [(100, 150), (300, 400)], hits)
        assert filtered[0].name == "A" and filtered[0].start == 120
        assert filtered[1].name == "B" and filtered[1].start == 300

    def test_load(self):
        matrices = load_matrices(PWM_PATH)
        assert matrices
        for matrix in matrices:
            assert isinstance(matrix, Matrix)


class TestFinder(unittest.TestCase):
    def setUp(self):
        self.region = DummyRegion()
        self.hit_record = DummyRecord(seq="AATTCCGGAAT", features=[self.region])
        self.region_hits = {
            1: [
                TFBSHit('A', 1, 'TestHit', 'TT', Confidence.WEAK, 1, 40.0, 41.7),
                TFBSHit('A', 1, 'TestHit', 'TT', Confidence.STRONG, 1, 40.0, 41.7),
                TFBSHit('B', 1, 'TestHit', 'TT', Confidence.MEDIUM, 1, 40.0, 41.7),
                TFBSHit('B', 2, 'TestHit', 'TT', Confidence.WEAK, 1, 19.0, 27.0),
            ],
        }
        self.results = TFBSFinderResults(record_id="dummy", pvalue=0.001, start_overlap=1,
                                         hits_by_region=self.region_hits)

    def test_filter_results(self):
        assert self.results.get_hits_by_region(2) == []
        assert self.results.get_hits_by_region(1) == sorted(self.region_hits[1])
        weak = [hit for hit in self.region_hits[1] if hit.confidence == Confidence.WEAK]
        assert self.results.get_hits_by_region(1, confidence=Confidence.WEAK) == weak
        other = [hit for hit in self.region_hits[1] if hit.confidence != Confidence.WEAK]
        assert self.results.get_hits_by_region(1, confidence=Confidence.MEDIUM, allow_better=True) == other
        assert sorted(other) == other

    def test_sequence_matches(self):
        input_seq = "TGTCCGATAGCAATCGACGT"
        consensus = "KKYYSSWWRRMMATCGNNNN"
        matches = get_sequence_matches(input_seq, consensus)
        assert matches == [True] * len(consensus)

        input_seq = "CGACTG"
        consensus = "KYSWRM"
        matches = get_sequence_matches(input_seq, consensus)
        assert matches == [False] * len(consensus)

    def test_invalid_matches(self):
        with self.assertRaisesRegex(ValueError, "must be the same length"):
            get_sequence_matches("A", "AAA")

    def test_generate_javascript_data(self):
        data = generate_javascript_data(self.hit_record, self.region, self.results)
        assert len(data) == 4

        assert data[0]['name'] == 'A'
        assert data[0]['start'] == 1
        assert data[0]['score'] == 40.0
        assert data[0]['confidence'] == 'strong'

        assert data[1]['name'] == 'B'
        assert data[1]['start'] == 1
        assert data[1]['score'] == 40.0
        assert data[1]['confidence'] == 'medium'

    def test_contained_gene(self):
        hit = {"start": 5, "end": 15}
        genes = [DummyCDS(start=0, end=3), DummyCDS(start=8, end=11), DummyCDS(start=13, end=18, strand=-1)]
        results = add_neighbouring_genes(hit, genes)
        assert results["left"] == {
            "location": 3,
            "name": genes[0].get_name(),
            "strand": 1,
        }
        assert results["mid"] == {
            "location": 8,
            "length": 3,
            "name": genes[1].get_name(),
            "strand": 1,
        }
        assert results["right"] == {
            "location": 13,
            "name": genes[2].get_name(),
            "strand": -1,
        }


class TestAreaFinding(unittest.TestCase):
    def setUp(self):
        self.region = DummyRegion(subregions=[DummySubRegion(start=100, end=500)])

    def test_leading(self):
        self.region.add_cds(DummyCDS(start=250, end=500))
        assert get_valid_areas(self.region, 0) == [(0, 250)]

    def test_trailing(self):
        self.region.add_cds(DummyCDS(start=0, end=250))
        assert get_valid_areas(self.region, 0) == [(250, 500)]

    def test_gap(self):
        self.region.add_cds(DummyCDS(start=0, end=200))
        self.region.add_cds(DummyCDS(start=300, end=500))
        assert get_valid_areas(self.region, 0) == [(200, 300)]

    def test_gaps_between_stops(self):
        self.region.add_cds(DummyCDS(start=0, end=180, strand=1))
        self.region.add_cds(DummyCDS(start=210, end=500, strand=-1))
        assert get_valid_areas(self.region, 0) == []

    def test_overlapping_ends(self):
        self.region.add_cds(DummyCDS(start=0, end=210, strand=1))
        self.region.add_cds(DummyCDS(start=200, end=500, strand=-1))
        assert get_valid_areas(self.region, 0) == []

    def test_overlapping_starts(self):
        self.region.add_cds(DummyCDS(start=0, end=210, strand=-1))
        self.region.add_cds(DummyCDS(start=200, end=500, strand=1))
        assert get_valid_areas(self.region, 50) == [(150, 260)]

    def test_overlapping_same_strand(self):
        for strand in [1, -1]:
            self.region._cdses.clear()
            self.region.add_cds(DummyCDS(start=0, end=210, strand=strand))
            self.region.add_cds(DummyCDS(start=200, end=500, strand=strand))
            assert get_valid_areas(self.region, 0) == []

    def test_overlap_size(self):
        for size in [0, 50, 100]:
            self.region._cdses.clear()
            self.region.add_cds(DummyCDS(start=50, end=200, strand=1))
            self.region.add_cds(DummyCDS(start=250, end=400, strand=-1))
            assert get_valid_areas(self.region, size) == [(0, 50 + size), (400 - size, 500)]
