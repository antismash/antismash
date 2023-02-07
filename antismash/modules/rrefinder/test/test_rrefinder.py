# License: GNU Affero General Public License v3 or later
# A copy of GNU AGPL v3 should have been included in this software package in LICENSE.txt.

# for test files, silence irrelevant and noisy pylint warnings
# pylint: disable=use-implicit-booleaness-not-comparison,protected-access,missing-docstring

from collections import defaultdict
import json as jsonlib
import unittest
from unittest.mock import patch

import antismash
from antismash.common.hmmer import HmmerResults
from antismash.common.secmet.features import FeatureLocation
from antismash.common.secmet.test.helpers import DummyRegion
from antismash.common.test.helpers import (
    DummyCandidateCluster,
    DummyCDS,
    DummyProtocluster,
    DummyRecord,
)
from antismash.common.test.test_hmmer import create_hmmer_hit as DummyHmmerHit
from antismash.config import build_config, destroy_config, get_config, update_config
from antismash.modules.rrefinder.html_output import will_handle
from antismash.modules.rrefinder.rrefinder import (
    RREFinderResults,
    check_hmm_hit,
    extract_rre_hits,
    filter_hits,
    gather_rre_candidates,
    run_rrefinder,
)
from antismash.modules.rrefinder.rre_domain import RREDomain, TOOL


class TestRREResults(unittest.TestCase):
    def setUp(self):
        self.hits = [
            DummyHmmerHit(locus_tag="a", protein_start=0, protein_end=50, score=38.0,
                          domain="RRE_type_A", identifier="RREFam001.1"),
            DummyHmmerHit(locus_tag="b", protein_start=0, protein_end=90, score=25.0,
                          domain="RRE_type_B", identifier="RREFam002.1")
        ]

        self.hits_by_cds = {hit.locus_tag: [hit] for hit in self.hits}
        self.hits_by_protocluster = {1: [self.hits[0].locus_tag],
                                     2: [self.hits[1].locus_tag]}
        self.tool = 'rrefinder'
        self.database = 'RREFam.hmm'
        self.detection = 'hmmscan'

        self.min_length = 50
        self.bitscore_cutoff = 25.0

        build_config([
            "--rre",
            "--rre-cutoff", f"{self.bitscore_cutoff:.1f}",
            "--rre-minlength", str(self.min_length),
        ], isolated=True, modules=antismash.get_all_modules())

        # Here a mock so that the call to its functions can be tracked
        self.record = unittest.mock.Mock()
        self.record.id = 'test_record'

    def tearDown(self):
        destroy_config()

    def create_results(self, record_id=None, cutoff=None, min_length=None,
                       hits_by_protocluster=None, hits_by_cds=None):
        record_id = record_id or self.record.id
        bitscore_cutoff = cutoff if cutoff is not None else self.bitscore_cutoff
        min_length = min_length or self.min_length
        hits_by_protocluster = hits_by_protocluster or self.hits_by_protocluster
        hits_by_cds = hits_by_cds or self.hits_by_cds

        return RREFinderResults(record_id, bitscore_cutoff, min_length,
                                hits_by_protocluster, hits_by_cds)

    def test_init(self):
        results = self.create_results()
        assert results.record_id == self.record.id
        assert results.bitscore_cutoff == self.bitscore_cutoff
        assert results.min_length == self.min_length
        assert results.hits_by_protocluster is self.hits_by_protocluster
        assert results.hits_by_cds is self.hits_by_cds
        assert results.tool == self.tool
        assert results.detection == self.detection
        assert results.database == self.database
        assert results.features

    def test_convert_hits_to_features(self):
        hit = DummyHmmerHit(location="[1000:2000]", locus_tag='a', identifier="RREFam123.1")

        result = self.create_results(hits_by_cds={hit.locus_tag: [hit]})
        assert len(result.features) == 1

        feature = result.features[0]
        assert isinstance(feature, RREDomain)
        assert feature.location == FeatureLocation(1000, 2000)
        assert feature.protein_location == FeatureLocation(hit.protein_start, hit.protein_end)
        assert feature.description == hit.description
        assert feature.domain == hit.domain
        assert feature.locus_tag == hit.locus_tag
        assert feature.domain_id == f"{self.tool}_{hit.locus_tag}_{hit.domain}.1"
        assert feature.database == self.database
        assert feature.detection == self.detection
        assert feature.score == hit.score
        assert feature.evalue == hit.evalue
        assert feature.label == hit.label
        assert feature.translation == hit.translation

    def test_json_conversion(self):
        results = self.create_results()

        json = jsonlib.loads(jsonlib.dumps(results.to_json()))

        regenerated = RREFinderResults.from_json(json, self.record)
        assert regenerated.hits_by_protocluster == results.hits_by_protocluster
        assert regenerated.hits_by_cds == results.hits_by_cds
        assert regenerated.bitscore_cutoff == results.bitscore_cutoff
        assert regenerated.min_length == results.min_length
        assert regenerated.record_id == results.record_id

    def test_from_json_wrong_record(self):
        json = self.create_results().to_json()
        json['record_id'] = 'another_record'
        assert not RREFinderResults.from_json(json, self.record)

    def test_from_json_wrong_schema_version(self):
        json = self.create_results().to_json()
        json['schema_version'] = 'not_a_valid_schema'
        assert not RREFinderResults.from_json(json, self.record)

    def test_from_json_lower_min_length(self):
        json = self.create_results().to_json()
        assert get_config().rre_min_length == 50
        update_config({"rre_min_length": 25})
        assert not RREFinderResults.from_json(json, self.record)

    def test_from_json_lower_bitscore(self):
        json = self.create_results().to_json()
        assert get_config().rre_cutoff == 25.
        update_config({"rre_cutoff": 15.})
        assert not RREFinderResults.from_json(json, self.record)

    def test_from_json_higher_min_length(self):
        json = self.create_results().to_json()
        assert get_config().rre_min_length == 50
        new = 80
        assert len(self.hits[0]) < new
        assert len(self.hits[1]) > new
        update_config({"rre_min_length": new})
        results = RREFinderResults.from_json(json, self.record)
        assert len(results.hits_by_cds) == 1
        assert results.hits_by_cds[self.hits[1].locus_tag] == [self.hits[1]]
        assert len(results.hits_by_protocluster) == 1
        assert results.hits_by_protocluster[2] == [self.hits[1].locus_tag]

    def test_from_json_higher_bitscore(self):
        json = self.create_results().to_json()
        assert get_config().rre_cutoff == 25.
        new = 35.
        assert self.hits[0].score > new
        assert self.hits[1].score < new
        update_config({"rre_cutoff": new})
        result = RREFinderResults.from_json(json, self.record)
        assert len(result.hits_by_cds) == 1
        assert result.hits_by_cds[self.hits[0].locus_tag] == [self.hits[0]]
        assert len(result.hits_by_protocluster) == 1
        assert result.hits_by_protocluster[1] == [self.hits[0].locus_tag]

    def test_add_to_incorrect_record(self):
        results = self.create_results(record_id=self.record.id)
        with self.assertRaisesRegex(ValueError, "Record to store in and record analysed don't match"):
            other = DummyRecord()
            other.id = self.record.id * 2
            results.add_to_record(other)

    def test_add_to_record(self):
        record = DummyRecord()
        results = self.create_results(record_id=record.id)

        assert not record.get_all_features()
        results.add_to_record(record)
        assert len(record.get_all_features()) == 2
        assert len(record.get_antismash_domains_by_tool(TOOL)) == 2


class TestRREFinder(unittest.TestCase):
    def setUp(self):
        self.hit_a1 = DummyHmmerHit(locus_tag='a', location='[1000:2000]',
                                    domain='RRE_type_A', evalue=0.1,
                                    protein_start=0, protein_end=50, score=38.0,
                                    identifier='RREFam001.1')
        self.hit_a2 = DummyHmmerHit(locus_tag='a', location='[1100:2100]',
                                    domain='RRE_type_A', evalue=0.1,
                                    protein_start=0, protein_end=40, score=38.0,
                                    identifier='RREFam001.1')
        self.hit_b1 = DummyHmmerHit(locus_tag='b', location='[3500:4500]',
                                    domain='RRE_type_B', evalue=0.2,
                                    protein_start=0, protein_end=90, score=25.0,
                                    identifier='RREFam002.1')
        self.hit_b2 = DummyHmmerHit(locus_tag='b', location='[3600:4600]',
                                    domain='RRE_type_B', evalue=0.2,
                                    protein_start=14, protein_end=50, score=42.0,
                                    identifier='RREFam003.1')
        self.hit_c = DummyHmmerHit(locus_tag='c', location='[200:400]',
                                   domain='RRE_type_D', evalue=0.2,
                                   protein_start=0, protein_end=120, score=6.0,
                                   identifier='RREFam004.1')

        self.candidates_per_protocluster = defaultdict(list)
        self.candidates_per_protocluster[1] = ['a', 'c']
        self.candidates_per_protocluster[2] = ['b']

        self.filtered_hits_by_protocluster = defaultdict(list)
        self.filtered_hits_by_protocluster[1] = ['a']
        self.filtered_hits_by_protocluster[2] = ['b']

        self.all_hits = [self.hit_a1, self.hit_a2, self.hit_b1, self.hit_b2, self.hit_c]
        self.fake_cds_info = dict((locus_tag, locus_tag) for locus_tag in ['a', 'b', 'c'])
        self.hits_by_cds = defaultdict(list)
        for hit in self.all_hits:
            self.hits_by_cds[hit.locus_tag].append(hit)
        self.filtered_hits_by_cds = {hit.locus_tag: [hit] for hit in [self.hit_a1, self.hit_b1]}

        self.bitscore_cutoff = 25.0
        self.max_evalue = 1
        self.min_length = 50
        self.record = self.make_dummy_record()
        self.record.id = 'test_record'
        self.non_ripp_record = self.make_nonripp_dummy_record()
        self.database = 'rre_database'
        self.tool = 'rrefinder'

        self.hmm_res = HmmerResults(self.record.id, self.max_evalue, self.bitscore_cutoff,
                                    self.database, self.tool, self.all_hits)

    def make_dummy_record(self):
        cds1 = DummyCDS(start=800, end=2150, locus_tag='a')
        cds2 = DummyCDS(start=3400, end=4700, locus_tag='b')
        cds3 = DummyCDS(start=150, end=450, locus_tag='c')

        p1 = DummyProtocluster(core_start=100, core_end=2200, neighbourhood_range=100,
                               product='lanthipeptide-class-i', product_category="RiPP")
        p2 = DummyProtocluster(core_start=3300, core_end=4800, neighbourhood_range=100,
                               product='thiopeptide', product_category="RiPP")

        dc1 = DummyCandidateCluster(clusters=[p1])
        dc2 = DummyCandidateCluster(clusters=[p2])

        region = DummyRegion(candidate_clusters=[dc1, dc2])
        return DummyRecord(seq='FAKESEQ'*1000, features=[cds1, cds2, cds3, p1, p2, region])

    def make_nonripp_dummy_record(self):
        region = DummyRegion()
        return DummyRecord(seq='FAKESEQ'*1000, features=[region])

    def test_will_handle(self):
        assert will_handle(["prod1", "prod2"], {"RiPP"})
        assert not will_handle(["prod1", "prod2"], {"nRiPP"})

    def test_check_hmm_hit(self):
        assert check_hmm_hit(self.hit_a1, self.min_length, self.bitscore_cutoff)
        assert not check_hmm_hit(self.hit_a2, self.min_length, self.bitscore_cutoff)
        assert check_hmm_hit(self.hit_b1, self.min_length, self.bitscore_cutoff)
        assert not check_hmm_hit(self.hit_b2, self.min_length, self.bitscore_cutoff)
        assert not check_hmm_hit(self.hit_c, self.min_length, self.bitscore_cutoff)

    def test_filter_hits(self):
        by_cds, by_proto = filter_hits(self.hits_by_cds, self.candidates_per_protocluster,
                                       self.min_length, self.bitscore_cutoff)
        assert self.filtered_hits_by_cds == by_cds
        assert self.filtered_hits_by_protocluster == by_proto

    def test_stricter_filters(self):
        by_cds, by_proto = filter_hits(self.hits_by_cds, self.candidates_per_protocluster, 1000, 10)
        assert not by_cds
        assert not by_proto

    def test_lenient_filters(self):
        by_cds, by_proto = filter_hits(self.hits_by_cds, self.candidates_per_protocluster, 0, 0)
        assert by_cds == self.hits_by_cds
        assert by_proto == self.candidates_per_protocluster

    def test_extract_rre_hits(self):
        assert self.hits_by_cds == extract_rre_hits(self.hmm_res)

    def test_gather_rre_candidates(self):
        candidates_per_protocluster, cds_info = gather_rre_candidates(self.record)
        assert sorted(self.candidates_per_protocluster) == sorted(candidates_per_protocluster)
        for key in self.fake_cds_info:
            assert key in cds_info
            assert key == cds_info[key].get_name()

    @patch('antismash.modules.rrefinder.rrefinder.run_hmmer')
    def test_run_rrefinder(self, mocked_function):
        mocked_function.return_value = self.hmm_res
        res_object = run_rrefinder(self.record, self.bitscore_cutoff, self.min_length, self.database)
        assert res_object.hits_by_protocluster == self.filtered_hits_by_protocluster
        assert res_object.hits_by_cds == self.filtered_hits_by_cds
        assert res_object.bitscore_cutoff == self.bitscore_cutoff
        assert res_object.min_length == self.min_length
        assert mocked_function.call_count == 1

    def test_run_rrefinder_no_ripps(self):
        res_object = run_rrefinder(self.non_ripp_record, self.bitscore_cutoff, self.min_length, self.database)
        assert res_object.hits_by_protocluster == {}
        assert res_object.hits_by_cds == {}

    @patch('antismash.modules.rrefinder.rrefinder.run_hmmer')
    def test_run_rrefinder_no_hits(self, mocked_function):
        mocked_function.return_value = HmmerResults(self.record.id, self.max_evalue, self.bitscore_cutoff,
                                                    self.database, self.tool, [])
        res_object = run_rrefinder(self.record, self.bitscore_cutoff, self.min_length, self.database)
        assert res_object.hits_by_protocluster == {}
        assert res_object.hits_by_cds == {}
