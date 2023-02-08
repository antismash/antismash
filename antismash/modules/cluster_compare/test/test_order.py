# License: GNU Affero General Public License v3 or later
# A copy of GNU AGPL v3 should have been included in this software package in LICENSE.txt.

# for test files, silence irrelevant and noisy pylint warnings
# pylint: disable=use-implicit-booleaness-not-comparison,protected-access,missing-docstring

import unittest

from antismash.common.test.helpers import DummyCDS
from antismash.modules.cluster_compare import ordering
from antismash.modules.cluster_compare.data_structures import Hit

from .test_data_structures import DummyReferenceCDS


def generate_hit(cds, ref_id):
    return Hit("rec", ref_id, cds, 50., 100., 50., 1e-8)


def score(segments, query_size, ref_size):
    result = ordering.score_segments(segments, query_size, ref_size)
    assert 0 <= result <= 1
    return result


class TestOrdering(unittest.TestCase):
    def setUp(self):
        self.refs = {}
        self.cdses = []
        for i in range(5):
            ref = DummyReferenceCDS(name=f"{i}", start=i*100)
            self.refs[ref.name] = ref
            self.cdses.append(DummyCDS(locus_tag=f"c{i}", start=i*100, end=i*109))

    def find(self, hits, cdses=None):
        return ordering.find_segments(hits, cdses or self.cdses, self.refs)

    def test_single(self):
        hits = {"0": generate_hit(self.cdses[0], "0")}
        assert self.find(hits) == [[(1, 0, True)]]

    def test_contiguous(self):
        hits = {
            "0": generate_hit(self.cdses[0], "0"),
            "1": generate_hit(self.cdses[1], "1"),
            "2": generate_hit(self.cdses[2], "2"),
        }
        assert self.find(hits) == [[(1, 0, True), (2, 1, True), (3, 2, True)]]

    def test_split(self):
        hits = {
            "0": generate_hit(self.cdses[0], "0"),
            "1": generate_hit(self.cdses[1], "1"),
            "4": generate_hit(self.cdses[2], "4"),
        }
        assert self.find(hits) == [[(1, 0, True), (2, 1, True)], [(3, 4, True)]]

    def test_scoring(self):
        hits = {
            "0": generate_hit(self.cdses[0], "0"),
            "1": generate_hit(self.cdses[1], "1"),
            "4": generate_hit(self.cdses[2], "4"),
        }
        segments = self.find(hits)
        expected = ordering.calculate_order_score(self.cdses, hits, self.refs)
        assert score(segments, len(self.cdses), len(self.refs)) == expected

    def test_empty(self):
        assert ordering.calculate_order_score([], {}, None) == 0.


class TestSegmentScoring(unittest.TestCase):
    def setUp(self):
        self.contiguous = [[(1, 2, True), (2, 3, True), (3, 4, True)]]
        self.single_split = [[(1, 2, True), (2, 3, True)], [(3, 10, True)]]
        self.multi_split = [[(1, 2, True)], [(2, 6, True)], [(3, 10, True)]]
        self.fully_reversed = [[(1, 10, False)], [(2, 6, False)], [(3, 2, False)]]

    def test_single(self):
        for segments in [[[(1, 2, True)]],
                         [[(1, 2, True), (2, 3, True)]],
                         ]:
            assert score(segments, len(segments[0]), len(segments[0])) == 1.

    def test_comparitive(self):
        assert score(self.contiguous, 3, 10) > score(self.single_split, 3, 10) > score(self.multi_split, 3, 10)
        assert score(self.contiguous, 5, 10) > score(self.single_split, 5, 10) > score(self.multi_split, 5, 10)

    def test_different_max(self):
        assert score(self.contiguous, 5, 10) < score(self.contiguous, 3, 10)

    def test_scaling(self):
        # numbers here are arbitrary, the point is to ensure that
        # a tiny continguous chunk of large number of genes is worth less than
        # a split chunk of a small number of genes
        tiny_contiguous = score([[(1, 2, True), (2, 3, True)]], 50, 50)
        intermediate = score([[(1, 2, True), (2, 3, True), (4, 5, True)]], 30, 30)
        split = score([[(1, 2, True), (3, 4, True)]], 4, 4)
        assert tiny_contiguous < intermediate
        assert intermediate < split

    def test_reversals(self):
        rev = self.multi_split[:]
        rev[1] = self.fully_reversed[1]
        assert score(rev, 5, 10) < score(self.multi_split, 5, 10)
        assert score(rev, 3, 10) < score(self.multi_split, 3, 10)
        # majority reversal should be accounted for
        rev = self.multi_split[:]
        rev[1] = self.multi_split[1]
        assert score(rev, 5, 10) == score(self.multi_split, 5, 10)
        assert score(self.fully_reversed, 5, 10) == score(self.multi_split, 5, 10)

    def test_gap_size(self):
        small = [[(1, 2, True)], [(2, 4, True)]]
        large = [[(1, 2, True)], [(2, 40, True)]]
        assert score(small, 2, 40) > score(large, 2, 40)
        assert score(small, 40, 40) > score(large, 40, 40)


class TestSegmentBuilding(unittest.TestCase):
    def find(self, pairings):
        result = ordering._build_segments_from_pairings(pairings)
        assert isinstance(result, list)
        if pairings:
            assert result
        return result

    def test_empty(self):
        assert self.find([]) == []

    def test_single_pair(self):
        assert self.find([(1, 5, True)]) == [[(1, 5, True)]]

    def test_bad_order(self):
        with self.assertRaisesRegex(ValueError, "pairs not ordered"):
            self.find([(2, 5, True), (1, 6, True)])

    def test_single_segment(self):
        pairs = [(1, 5, True), (2, 6, True)]
        result = self.find(pairs)
        assert result == [pairs]

    def test_reverse_strand_single_segment(self):
        pairs = [(1, 5, False), (2, 4, False)]
        result = self.find(pairs)
        assert result == [pairs]

    def test_split_segments(self):
        pairs = [(1, 5, True), (3, 6, True)]
        result = self.find(pairs)
        assert result == [[pairs[0]], [pairs[1]]]

    def test_split_reverse_segments(self):
        pairs = [(1, 5, False), (3, 6, False)]
        result = self.find(pairs)
        assert result == [[pairs[0]], [pairs[1]]]

    def test_split_strand(self):
        pairs = [(1, 5, True), (2, 6, False), (3, 7, False)]
        result = self.find(pairs)
        assert result == [[pairs[0]], pairs[1:]]
