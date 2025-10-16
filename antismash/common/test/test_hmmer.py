# License: GNU Affero General Public License v3 or later
# A copy of GNU AGPL v3 should have been included in this software package in LICENSE.txt.

# for test files, silence irrelevant and noisy pylint warnings
# pylint: disable=use-implicit-booleaness-not-comparison,protected-access,missing-docstring

from collections import defaultdict
from dataclasses import FrozenInstanceError, dataclass
import os
from tempfile import NamedTemporaryFile
import unittest
from unittest.mock import patch

from antismash.common import hmmer, json, subprocessing
from antismash.common.hmmer import (
    ensure_database_pressed,
    HmmerHit,
    HmmerResults,
    remove_overlapping,
)

from .helpers import DummyRecord


def create_hmmer_hit(location="[500:700]", label="ref_name", locus_tag="locus",
                     domain="hit", evalue=1e-10, score=40.5, translation=None,
                     identifier="TESTID_0005", description="some description",
                     protein_start=9, protein_end=None):
    if translation is None:
        if protein_end is None:
            protein_end = 13
        translation = "M" * (protein_end - protein_start)

    protein_end = protein_end or (protein_start + len(translation))

    return HmmerHit(
        location=location,
        label=label,
        locus_tag=locus_tag,
        domain=domain,
        evalue=evalue,
        score=score,
        translation=translation,
        identifier=identifier,
        description=description,
        protein_start=protein_start,
        protein_end=protein_end
    )


class TestHmmerHit(unittest.TestCase):
    def test_immutable(self):
        hit = create_hmmer_hit()
        with self.assertRaisesRegex(FrozenInstanceError, "cannot assign to field 'invalid'"):
            hit.invalid = "test"
        assert hit.protein_start
        with self.assertRaisesRegex(FrozenInstanceError, "cannot assign to field 'protein_start'"):
            hit.protein_start = 7

    def test_length(self):
        hit = create_hmmer_hit(protein_start=9, protein_end=13)
        assert len(hit) == 4
        hit = create_hmmer_hit(protein_start=1, protein_end=23)
        assert len(hit) == 22

    def test_invalid_positions(self):
        with self.assertRaisesRegex(ValueError, "inverted start and end"):
            create_hmmer_hit(protein_start=9, protein_end=9)
        with self.assertRaisesRegex(ValueError, "inverted start and end"):
            create_hmmer_hit(protein_start=5, protein_end=3)

    def test_mismatching_lengths(self):
        hit = create_hmmer_hit(protein_start=5, protein_end=15, translation="A"*10)
        assert len(hit) == len(hit.translation) == hit.protein_end - hit.protein_start
        with self.assertRaisesRegex(ValueError, "translation length does not match"):
            create_hmmer_hit(protein_start=5, protein_end=15, translation="A"*4)

    def test_json_conversion(self):
        hit = create_hmmer_hit()
        as_json = json.loads(json.dumps(hit.to_json()))
        assert hit == HmmerHit.from_json(as_json)


class TestOverlaps(unittest.TestCase):
    @dataclass(frozen=True)
    class Dummy:
        identifier: str
        protein_start: int
        protein_end: int
        score: float = 10.

        def __repr__(self):
            # useful when tests fail
            return self.identifier  # pragma: no cover

        def __len__(self):
            return self.protein_end - self.protein_start

    def setUp(self):
        self.cutoffs = defaultdict(lambda: 5.)
        # these names match the effective scores, not the cutoff themselves
        self.cutoffs.update({"low": 10., "med": 5., "high": 1.})

    def test_no_overlap(self):
        hit_a = self.Dummy("A", 50, 100)
        hit_b = self.Dummy("B", 150, 200)
        assert remove_overlapping([hit_a, hit_b], self.cutoffs) == [hit_a, hit_b]
        assert remove_overlapping([hit_b, hit_a], self.cutoffs) == [hit_a, hit_b]

    def test_limit(self):
        hit_a = self.Dummy("A", 50, 200)
        hit_b = self.Dummy("B", 150, 200)
        assert len(remove_overlapping([hit_a, hit_b], self.cutoffs, overlap_limit=100)) == 2
        assert len(remove_overlapping([hit_a, hit_b], self.cutoffs, overlap_limit=10)) == 1

    def test_single_overlap(self):
        hit_a = self.Dummy("low", 50, 150)
        hit_b = self.Dummy("high", 100, 200)
        assert remove_overlapping([hit_a, hit_b], self.cutoffs, overlap_limit=10) == [hit_b]
        assert remove_overlapping([hit_b, hit_a], self.cutoffs, overlap_limit=10) == [hit_b]

    def test_multiple_rising(self):
        low = self.Dummy("low", 50, 150)
        med = self.Dummy("med", 120, 220)
        high = self.Dummy("high", 200, 300)
        assert remove_overlapping([low, med, high], self.cutoffs, overlap_limit=10) == [low, high]

    def test_multiple_falling(self):
        high = self.Dummy("high", 50, 150)
        med = self.Dummy("med", 120, 220)
        low = self.Dummy("low", 200, 300)
        assert remove_overlapping([low, med, high], self.cutoffs, overlap_limit=10) == [high, low]

    def test_troughs(self):
        first = self.Dummy("high", 50, 150)
        mid = self.Dummy("med", 120, 220)
        last = self.Dummy("high", 200, 300)
        assert remove_overlapping([first, mid, last], self.cutoffs, overlap_limit=10) == [first, last]

    def test_peaks(self):
        first = self.Dummy("low", 50, 150)
        mid = self.Dummy("med", 120, 220)
        last = self.Dummy("low", 200, 300)
        assert remove_overlapping([first, mid, last], self.cutoffs, overlap_limit=10) == [mid]

    def test_tiebreaker_length(self):
        hit_a = self.Dummy("a", 50, 125)
        hit_b = self.Dummy("b", 100, 200)
        # should take the longest in case of hit_a tie
        assert remove_overlapping([hit_a, hit_b], self.cutoffs, overlap_limit=10) == [hit_b]
        hit_a = self.Dummy("a", 50, 175)
        assert remove_overlapping([hit_a, hit_b], self.cutoffs, overlap_limit=10) == [hit_a]

    def test_tiebreaker_start(self):
        a = self.Dummy("a", 50, 150)
        b = self.Dummy("b", 100, 200)
        # should take the longest in case of a tie
        assert remove_overlapping([a, b], self.cutoffs, overlap_limit=10) == [a]
        assert remove_overlapping([b, a], self.cutoffs, overlap_limit=10) == [a]


class TestResults(unittest.TestCase):
    def setUp(self):
        self.low_score = create_hmmer_hit(score=50., evalue=1e-20)
        self.high_score = create_hmmer_hit(score=75., evalue=1e-20)
        self.low_evalue = create_hmmer_hit(score=60., evalue=1e-40)
        self.high_evalue = create_hmmer_hit(score=60., evalue=1e-10)
        hits = [self.low_score, self.high_score, self.low_evalue, self.high_evalue]

        self.results = HmmerResults("name", 1e-1, 25., "db", "tool", hits)

    def test_refilter_higher_score(self):
        assert len(self.results.hits) == 4
        self.results.refilter(1e-10, 55.)
        assert len(self.results.hits) == 3
        self.results.refilter(1e-10, 65.)
        assert len(self.results.hits) == 1
        self.results.refilter(1e-10, 85.)
        assert len(self.results.hits) == 0

    def test_refilter_lower_evalue(self):
        assert len(self.results.hits) == 4
        self.results.refilter(1e-15, 50.)
        assert len(self.results.hits) == 3
        self.results.refilter(1e-25, 50.)
        assert len(self.results.hits) == 1
        self.results.refilter(1e-50, 50.)
        assert len(self.results.hits) == 0

    def test_refilter_lower_score(self):
        with self.assertRaisesRegex(ValueError, "more lenient"):
            self.results.refilter(self.results.evalue, self.results.score / 2)

    def test_refilter_higher_evalue(self):
        with self.assertRaisesRegex(ValueError, "more lenient"):
            self.results.refilter(self.results.evalue * 2, self.results.score)

    def test_json_conversion(self):
        as_string = json.dumps(self.results.to_json())
        regenerated = HmmerResults.from_json(json.loads(as_string), DummyRecord(record_id="name"))
        assert isinstance(regenerated, HmmerResults)
        assert regenerated.hits == self.results.hits
        assert regenerated.evalue == self.results.evalue
        assert regenerated.score == self.results.score
        assert regenerated.database == self.results.database
        assert regenerated.tool == self.results.tool


class TestProfileAggregation(unittest.TestCase):
    def setUp(self):
        self.aggregate = NamedTemporaryFile()
        self.first = NamedTemporaryFile()
        self.second = NamedTemporaryFile()
        for temp, content in [(temp, content) for temp, content in zip([self.first, self.second], "AB")]:
            temp.write(content.encode())
            temp.flush()

    def tearDown(self):
        self.aggregate.close()
        self.first.close()
        self.second.close()

    def test_aggregation(self):
        # aggregate is oldest
        with patch.object(os.path, "getmtime", side_effect=[1, 2, 3]):
            with patch.object(hmmer, "ensure_database_pressed", return_value=[]):
                errors = hmmer.aggregate_profiles(self.aggregate.name, [self.first.name, self.second.name])
        assert not errors
        # content should have been written
        with open(self.aggregate.name, "r") as handle:
            assert handle.read() == "AB"

        self.aggregate.close()
        self.aggregate = NamedTemporaryFile()

        # aggregate is newest
        with patch.object(os.path, "getmtime", side_effect=[3, 1, 2]):
            with patch.object(hmmer, "ensure_database_pressed", return_value=[]):
                errors = hmmer.aggregate_profiles(self.aggregate.name, [self.first.name, self.second.name])
        assert not errors
        # nothing should have been written
        assert not self.aggregate.read()


class TestPressed(unittest.TestCase):
    def test_caught(self):
        with self.assertRaisesRegex(ValueError, "set to mount at runtime"):
            ensure_database_pressed("/path/mounted_at_runtime/pfam")


class TestRunner(unittest.TestCase):
    def test_empty_features(self):
        record = DummyRecord()
        with patch.object(subprocessing.hmmscan, "run_hmmscan") as patched:
            with patch.object(os.path, "exists", return_value=True):
                results = hmmer.run_hmmer(record, features=[], max_evalue=1.,
                                          min_score=0, database="", tool="dummy")
            patched.assert_not_called()
        assert not results.hits
