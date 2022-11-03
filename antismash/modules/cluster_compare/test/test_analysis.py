# License: GNU Affero General Public License v3 or later
# A copy of GNU AGPL v3 should have been included in this software package in LICENSE.txt.

# for test files, silence irrelevant and noisy pylint warnings
# pylint: disable=use-implicit-booleaness-not-comparison,protected-access,missing-docstring

import math
import unittest
from unittest.mock import Mock

from antismash.common.secmet.features.cdscollection import CDSCollection, FeatureLocation
from antismash.common.test.helpers import DummyCDS
from antismash.modules.cluster_compare.analysis import (
    blast_parse,
    calculate_identity_score as ident,
    calculate_protocluster_ranking,
    convert_to_references as convert,
    filter_by_query_area,
    filter_by_reference_protocluster,
    Hit,
    trim_to_best_hit as trim,
    parse_hit,
)
from antismash.modules.cluster_compare.data_structures import (
    ReferenceCDS,
    ReferenceProtocluster,
    ReferenceRecord,
    ReferenceRegion,
)


class DummyHit(Hit):
    def __init__(self, pid=.5, cds=None, ref_name="ref_cds_name"):
        super().__init__("ref_rec", ref_name, cds or DummyCDS(), pid, 1234., 105., 1e-8)


class TestCalculations(unittest.TestCase):
    def test_identity_extremes(self):
        # this might not seem necessary, but actual calculation could result in
        # values outside the expected range, so ensure they're caught
        assert ident([DummyHit(pid=0.0)], 1) == 0.
        assert ident([DummyHit(pid=1e-9)], 1) == 0.
        assert ident([DummyHit(pid=100.0 - 1e-9)], 1) == 1.
        assert ident([DummyHit(pid=100.0)], 1) == 1.

    def test_identity_averaging(self):
        # the input to the core scoring function should use the average PID
        hit = DummyHit(pid=80.)
        for i in range(1, 4):
            self.assertAlmostEqual(ident([hit]*i, i), 0.615, delta=1e-3)
        self.assertAlmostEqual(ident([DummyHit(pid=90.), DummyHit(pid=70.)], 2), 0.615, delta=1e-3)

    def test_identity_no_hits(self):
        assert ident(iter([]), 0) == 0.
        # even if we lie about the number of hits
        assert ident(iter([]), 5) == 0

    def test_ranking(self):
        def rank_calc(score):
            scorer = Mock()
            scorer.final_score = score
            return calculate_protocluster_ranking(scorer)

        assert rank_calc(1.) == 1.
        assert rank_calc(0.) == 0.

        assert rank_calc(.1) < rank_calc(.4) < rank_calc(.8)


class TestTrimToBest(unittest.TestCase):
    def setUp(self):
        cds_a = DummyCDS(locus_tag="A")
        cds_b = DummyCDS(locus_tag="B")
        self.a_high = DummyHit(pid=80, cds=cds_a)
        self.a_med = DummyHit(pid=50, cds=cds_a)
        self.a_low = DummyHit(pid=30, cds=cds_a)
        self.b_high = DummyHit(pid=80, cds=cds_b)
        self.b_low = DummyHit(pid=40, cds=cds_b)

    def test_single_set(self):
        assert trim({"A": [self.b_low, self.a_high, self.a_low]}) == {"A": self.a_high}

    def test_removed_completely(self):
        assert trim({"A": [self.a_high], "B": [self.a_low]}) == {"A": self.a_high}
        assert trim({"A": [self.b_low], "B": [self.b_high]}) == {"B": self.b_high}

    def test_runner_up(self):
        hits = {"A": [self.a_med, self.b_high], "B": [self.b_high, self.a_low]}
        assert trim(hits) == {"A": self.a_med, "B": self.b_high}


class TestParsing(unittest.TestCase):
    def test_parse_line(self):
        cds = DummyCDS(start=1, end=500)
        line = "1\trefrec1.1|21\t53.8\t234\t97\t7\t14\t240\t2\t231\t1.3e-55\t214.5"
        hit = parse_hit(line, {1: cds})
        assert isinstance(hit, Hit)
        assert hit.cds is cds
        assert hit.reference_record == "refrec1.1"
        assert hit.reference_id == "21"
        assert hit.percent_identity == 53.8
        assert hit.evalue == 1.3e-55
        assert hit.blast_score == 214.5
        assert math.floor(hit.percent_coverage) == 46

    def test_parse_all(self):
        blast_output = ("1\trefrec1.1|21\t53.8\t234\t97\t7\t14\t240\t2\t231\t1.3e-55\t214.5\n"
                        "2\trefrec2|3\t87.7\t130\t16\t0\t1\t130\t1\t130\t2.0e-66\t249.6\n")
        cdses = {1: DummyCDS(start=1, end=500), 2: DummyCDS(start=600, end=700)}
        by_cds, by_ref = blast_parse(blast_output, cdses, min_seq_coverage=0, min_perc_identity=0)
        assert len(by_cds) == 2
        assert len(by_ref) == 2

        hit = by_ref["refrec1.1"]["21"][0]
        assert hit.percent_identity == 53.8
        assert hit.cds == cdses[1]
        assert hit is by_cds[cdses[1]]["refrec1.1"][0]

        hit = by_ref["refrec2"]["3"][0]
        assert hit.percent_identity == 87.7
        assert hit.cds == cdses[2]
        assert hit is by_cds[cdses[2]]["refrec2"][0]

        by_cds, by_ref = blast_parse(blast_output, cdses, min_seq_coverage=0, min_perc_identity=80)
        assert len(by_ref) == len(by_cds) == 1
        by_cds, by_ref = blast_parse(blast_output, cdses, min_seq_coverage=50, min_perc_identity=0)
        assert len(by_ref) == len(by_cds) == 1


class TestFiltering(unittest.TestCase):
    def test_by_query(self):
        area = CDSCollection(FeatureLocation(100, 500, 1), "test")
        contained_cds = DummyCDS(start=150, end=250)
        area.add_cds(contained_cds)
        missed_cds = DummyCDS()
        ref = ReferenceRecord("test", [], {})
        hits = {ref: {"a": [DummyHit(cds=contained_cds), DummyHit(cds=missed_cds)]}}
        trimmed = filter_by_query_area(area, hits)
        assert trimmed == {ref: {"a": [hits[ref]["a"][0]]}}

    def test_by_reference(self):
        ref_record = ReferenceRecord("test", [], {})
        ref_cds = ReferenceCDS("a", "other", {}, FeatureLocation(50, 450))
        ref_proto = ReferenceProtocluster("test", 1, 500, {}, {"a": ref_cds}, [], "product")
        hits = {ref_record: {"a": [DummyHit()], "b": [DummyHit()]}}
        trimmed = filter_by_reference_protocluster(ref_proto, hits)
        assert "b" not in trimmed[ref_record]


class TestConversion(unittest.TestCase):
    def test_convert(self):
        loc = FeatureLocation(5, 50, 1)
        records = {}
        hits = {}
        for i in "ab":
            cds = ReferenceCDS(f"CDS{i}", "other", {}, loc)
            region = ReferenceRegion(i, loc.start, loc.end, [], {cds.name: cds}, [], {"1": cds.name}, "desc", "org")
            records[i] = ReferenceRecord(i, [region], {"1": cds.name})
            hits[i] = {"1": [DummyHit(ref_name=cds.name)]}
        result = convert(hits, records)
        assert list(result) == [record.regions[0] for record in records.values()]
        assert result[records["a"].regions[0]] == {"CDSa": hits["a"]["1"]}
