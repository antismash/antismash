# License: GNU Affero General Public License v3 or later
# A copy of GNU AGPL v3 should have been included in this software package in LICENSE.txt.

# for test files, silence irrelevant and noisy pylint warnings
# pylint: disable=use-implicit-booleaness-not-comparison,protected-access,missing-docstring

import unittest
import warnings

from antismash.common import json, path
from antismash.common.test.helpers import DummyHMMResult
import antismash.common.hmmscan_refinement as refinement

# Don't display the SearchIO experimental warning, we know this.
with warnings.catch_warnings():
    warnings.simplefilter("ignore")
    from Bio import SearchIO


class TestHMMResult(unittest.TestCase):
    def test_construction(self):
        strings = refinement.HMMResult("dummy_hit", "1", "5", "3E-10", "53.5")
        correct = refinement.HMMResult("dummy_hit", 1, 5, 3e-10, 53.5)
        for result in [strings, correct]:
            assert result.hit_id == "dummy_hit"
            assert result.query_start == 1
            assert result.query_end == 5
            assert result.evalue == 3e-10
            assert result.bitscore == 53.5

    def test_length(self):
        result = refinement.HMMResult("dummy_hit", 1, 5, 3e-10, 53.5, 1, 6)
        assert result.query_length == 4
        assert result.hit_length == 5

    def test_merge(self):
        first = refinement.HMMResult("dummy_hit", 1, 5, 3e-10, 53.5, 1, 6)
        second = refinement.HMMResult("dummy_hit", 15, 25, 3e-20, 73.5, 20, 30)
        for merged in [first.merge(second), second.merge(first)]:
            assert merged.hit_id == "dummy_hit"
            assert merged.hit_start == 1
            assert merged.hit_end == 30
            assert merged.query_start == 1
            assert merged.query_end == 25
            assert merged.evalue == 3e-20
            assert merged.bitscore == 73.5

    def test_json_conversion(self):
        result = refinement.HMMResult("dummy_hit", 1, 5, 3e-10, 53.5)
        data = result.to_json()
        assert data == {'bitscore': 53.5,
                        'evalue': 3e-10,
                        'hit_id': 'dummy_hit',
                        'hit_start': 1,
                        'hit_end': 1,
                        'query_end': 5,
                        'query_start': 1}
        regenerated = refinement.HMMResult.from_json(data)
        assert regenerated.hit_id == "dummy_hit"
        assert regenerated.query_start == 1
        assert regenerated.query_end == 5
        assert regenerated.evalue == 3e-10
        assert regenerated.bitscore == 53.5

    def test_str_conversion(self):
        result = refinement.HMMResult("dummy_hit", 1, 5, 3e-10, 53.5, 1, 6)
        assert str(result) == "HMMResult(dummy_hit, hit_start=1, hit_end=6, " \
                              "query_start=1, query_end=5, evalue=3e-10, bitscore=53.5)"
        outer = refinement.HMMResult("other", 1, 5, 3e-10, 53.5, internal_hits=[result])
        assert str(outer).endswith(", subtypes=[dummy_hit])")

    def test_equality(self):
        first = refinement.HMMResult("dummy_hit", 1, 5, 3e-10, 53.5)
        second = refinement.HMMResult("dummy_hit", 1, 5, 3e-10, 53.5)
        assert first == second and first is not second
        second._hit_id = "dummy"
        assert first != second
        second._hit_id = first._hit_id
        second._evalue /= 10
        assert first != second

    def test_query_containment(self):
        outer = DummyHMMResult(start=5, end=20)
        # check self containment
        assert outer.query_is_contained_by(outer)
        # check right bounds
        inner = DummyHMMResult(start=10, end=20)
        assert inner.query_is_contained_by(outer)
        assert not outer.query_is_contained_by(inner)
        # check left bounds
        inner = DummyHMMResult(start=5, end=15)
        assert inner.query_is_contained_by(outer)
        assert not outer.query_is_contained_by(inner)
        # check completely within
        inner = DummyHMMResult(start=10, end=15)
        assert inner.query_is_contained_by(outer)
        assert not outer.query_is_contained_by(inner)

    def test_hit_containment(self):
        outer = DummyHMMResult(hit_start=5, hit_end=20)
        # check self containment
        assert outer.hit_is_contained_by(outer)
        # check right bounds
        inner = DummyHMMResult(hit_start=10, hit_end=20)
        assert inner.hit_is_contained_by(outer)
        assert not outer.hit_is_contained_by(inner)
        # check left bounds
        inner = DummyHMMResult(hit_start=5, hit_end=15)
        assert inner.hit_is_contained_by(outer)
        assert not outer.hit_is_contained_by(inner)
        # check completely within
        inner = DummyHMMResult(hit_start=10, hit_end=15)
        assert inner.hit_is_contained_by(outer)
        assert not outer.hit_is_contained_by(inner)

    def test_query_overlaps(self):
        center = DummyHMMResult(start=10, end=20)
        # Neigbouring to the left
        floating = DummyHMMResult(start=1, end=10)
        assert not center.query_overlaps_with(floating)
        # Overlapping the left bound
        floating = DummyHMMResult(start=1, end=15)
        assert center.query_overlaps_with(floating)
        assert center.query_overlaps_with(floating, min_overlap=5)
        assert not center.query_overlaps_with(floating, min_overlap=6)
        # Overlapping the right bound
        floating = DummyHMMResult(start=15, end=25)
        assert center.query_overlaps_with(floating)
        assert center.query_overlaps_with(floating, min_overlap=5)
        assert not center.query_overlaps_with(floating, min_overlap=6)
        # Neigbouring to the right
        floating = DummyHMMResult(start=20, end=30)
        assert not center.query_overlaps_with(floating)

    def test_hit_overlaps(self):
        center = DummyHMMResult(hit_start=10, hit_end=20)
        # Neigbouring to the left
        floating = DummyHMMResult(hit_start=1, hit_end=10)
        assert not center.hit_overlaps_with(floating)
        # Overlapping the left bound
        floating = DummyHMMResult(hit_start=1, hit_end=15)
        assert center.hit_overlaps_with(floating)
        assert center.hit_overlaps_with(floating, min_overlap=5)
        assert not center.hit_overlaps_with(floating, min_overlap=6)
        # Overlapping the right bound
        floating = DummyHMMResult(hit_start=15, hit_end=25)
        assert center.hit_overlaps_with(floating)
        assert center.hit_overlaps_with(floating, min_overlap=5)
        assert not center.hit_overlaps_with(floating, min_overlap=6)
        # Neigbouring to the right
        floating = DummyHMMResult(hit_start=20, hit_end=30)
        assert not center.hit_overlaps_with(floating)

    def test_hashability(self):
        first = refinement.HMMResult("dummy_hit", 1, 5, 3e-10, 53.5)
        second = refinement.HMMResult("dummy_hit", 1, 5, 3e-10, 53.5)
        different = refinement.HMMResult("dummy_hit", 1, 5, 3e-1, 53.5)

        assert hash(first) == hash(second) and first is not second
        assert hash(first) != hash(different)

        used = {}
        used[first] = 1
        used[second] = 2
        used[different] = 3

        assert used == {first: 2, different: 3}

    def test_nesting(self):
        inner = refinement.HMMResult("in", 5, 10, 3e-10, 53.5)
        assert inner.internal_hits == tuple()
        mid = refinement.HMMResult("mid", 3, 12, 3e-10, 53.5, internal_hits=[inner])
        assert mid.internal_hits == (inner,)
        outer = refinement.HMMResult("out", 3, 15, 3e-10, 53.5, internal_hits=[mid])
        assert outer.internal_hits == (mid,)
        assert outer.internal_hits[0].internal_hits[0] is inner

        reconstructed = refinement.HMMResult.from_json(json.loads(json.dumps(outer.to_json())))
        assert reconstructed == outer
        assert reconstructed.internal_hits == outer.internal_hits
        assert reconstructed.internal_hits[0].internal_hits == mid.internal_hits

        # and names check
        assert outer.detailed_names == [outer.hit_id, mid.hit_id, inner.hit_id]
        outer.add_internal_hits([inner])
        assert outer.detailed_names == [outer.hit_id]

    def test_nesting_additions(self):
        inner = refinement.HMMResult("in", 5, 7, 3e-10, 53.5)
        other = refinement.HMMResult("in", 8, 10, 3e-10, 53.5)
        outer = refinement.HMMResult("mid", 3, 12, 3e-10, 53.5, internal_hits=[inner])
        outer.add_internal_hits([other])
        assert outer.internal_hits == (inner, other)


class TestRefinement(unittest.TestCase):
    def setUp(self):
        self.gene_id = "SCO1217"
        datafile = path.get_full_path(__file__, 'data', 'hmmscan', 'SCO1217.hmmer3.txt')
        self.results = list(SearchIO.parse(datafile, 'hmmer3-text'))
        self.hmm_lengths = {"SMCOG1003:sensor_histidine_kinase": 570,
                            "SMCOG1048:sensor_histidine_kinase": 417,
                            "SMCOG1237:transposase": 387}
        self.hit_ids = set(self.hmm_lengths)

    # TODO: create a test with multiple genes hit

    def test_gather(self):
        # ensure the test will function as expected
        for result in self.results:
            for hsps in result.hsps:
                assert hsps.query_id == self.gene_id

        gathered = refinement.gather_by_query(self.results)
        # make sure we only found the one we're interested
        assert len(gathered) == 1
        assert [self.gene_id] == list(gathered)

        for result in gathered[self.gene_id]:
            assert isinstance(result, refinement.HMMResult)

    def test_merges(self):
        results = refinement.gather_by_query(self.results)[self.gene_id]
        results = sorted(list(results), key=lambda result: result.query_start)
        hit_ids = set(result.hit_id for result in results)
        assert len(hit_ids) == 3
        assert hit_ids == self.hit_ids
        # 2 hits each for 1048 and 1237, for 5 hits total
        assert len(results) == 5
        new = refinement._merge_domain_list(results, self.hmm_lengths)
        # after merging, only one hit should remain for each
        assert len(new) == 3
        # and the original list should be untouched
        assert len(results) == 5

    def test_incomplete_removal(self):
        results = refinement.gather_by_query(self.results)[self.gene_id]
        results = sorted(list(results), key=lambda result: result.query_start)
        assert len(results) == 5
        # ensure they're all too short to be caught
        for result in results:
            assert result.query_length / self.hmm_lengths[result.hit_id] < 1
        new = refinement.remove_incomplete(results, self.hmm_lengths)
        # ensure all were removed
        assert not new
        # and original list untouched
        assert len(results) == 5

        longest = 0
        for result in results:
            proportional_length = result.query_length / self.hmm_lengths[result.hit_id]
            if proportional_length > longest:
                longest = proportional_length

        assert longest < 1./3.
        # ensure the fallback works as intended
        new = refinement.remove_incomplete(results, self.hmm_lengths)
        # ensure all were removed
        assert not new

        new = refinement.remove_incomplete(results, self.hmm_lengths, fallback=longest - 0.01)
        # ensure the longest, and longer than the fallback, remain
        assert len(new) == 1
        assert new[0].query_length / self.hmm_lengths[new[0].hit_id] == longest

        # change the fallback to 0 and ensure only one comes back
        new = refinement.remove_incomplete(results, self.hmm_lengths, fallback=0.)
        assert len(new) == 1
        assert new[0].query_length / self.hmm_lengths[new[0].hit_id] == longest

    def test_incomplete_regulator(self):
        results = refinement.gather_by_query(self.results)[self.gene_id]
        results = sorted(list(results), key=lambda result: result.query_start)
        assert len(results) == 5
        regulator_id = "DUMMY:some_regulator_desc"
        regulator_result = refinement.HMMResult(regulator_id, 1, 2, 1e-10, 1)
        results.append(regulator_result)
        new_lengths = dict(self.hmm_lengths)
        new_lengths[regulator_id] = regulator_result.query_length * 100  # always big
        # set the thresholds to be unreachable
        new = refinement.remove_incomplete(results, new_lengths, threshold=2., fallback=2.)
        # ensure the tiny, but present, regulator is still in the list
        assert len(new) == 1
        assert new[0].hit_id == regulator_id

    def test_overlaps(self):
        first = refinement.HMMResult("dummy_hit", 1, 10, 3e-10, 53.5)
        second = refinement.HMMResult("dummy_hit", 10, 20, 3e-20, 73.5)
        results = [first, second]
        assert len(refinement._remove_overlapping(results, {"dummy_hit": 20})) == 2
        assert len(refinement._remove_overlapping(results, {"dummy_hit": 100})) == 2
        first._query_end = 16
        assert len(refinement._remove_overlapping(results, {"dummy_hit": 20})) == 1
        assert len(refinement._remove_overlapping(results, {"dummy_hit": 100})) == 2

        first._query_end = 13
        assert len(refinement._remove_overlapping(results, {"dummy_hit": 10})) == 1
        assert len(refinement._remove_overlapping(results, {"dummy_hit": 100})) == 2

    def test_combined(self):
        results = refinement.refine_hmmscan_results(self.results, self.hmm_lengths)
        assert len(results) == 1
        assert len(results[self.gene_id]) == 1
        best = results[self.gene_id][0]
        assert best.hit_id == "SMCOG1048:sensor_histidine_kinase"
        assert best.evalue == 3.6e-13
        assert best.bitscore == 43.5
        assert best.query_start == 91
        assert best.query_end == 390
