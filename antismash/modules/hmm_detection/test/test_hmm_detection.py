# License: GNU Affero General Public License v3 or later
# A copy of GNU AGPL v3 should have been included in this software package in LICENSE.txt.

# for test files, silence irrelevant and noisy pylint warnings
# pylint: disable=no-self-use,protected-access,missing-docstring

import os
import unittest

from Bio.Seq import Seq
from minimock import mock, restore

from antismash.common.deprecated import FeatureLocation
from antismash.common.test.helpers import DummyRecord, DummyCDS, DummyFeature, DummyCluster
from antismash.common.secmet import Record, CDSFeature, Feature
import antismash.common.path as path
from antismash.config import get_config, update_config, destroy_config, args
import antismash.modules.hmm_detection as core
from antismash.modules.hmm_detection import hmm_detection, rule_parser, signatures


class FakeHSP(object):  # pylint: disable=too-few-public-methods
    "class for generating a HSP like datastructure"
    def __init__(self, query_id, hit_id, hit_start, hit_end, bitscore, evalue):
        self.query_id = query_id
        self.hit_id = hit_id
        self.hit_start = hit_start
        self.hit_end = hit_end
        self.bitscore = bitscore
        self.evalue = evalue

    def __repr__(self):
        return "FakeHSP({})".format(str(vars(self)))


class TestArgs(unittest.TestCase):
    def setUp(self):
        self.parser = args.build_parser(modules=[core])

    def tearDown(self):
        destroy_config()

    def build_options(self, options):
        destroy_config()
        options = self.parser.parse_args(options)
        return update_config(options)

    def test_args(self):
        # argparse raises SystemExit, because it's missing the following option
        with self.assertRaises(SystemExit):
            self.build_options(["--enable"])
        options = self.build_options(["--enable", "t1pks,other"])
        assert options.enabled_cluster_types == ["t1pks", "other"]
        options = self.build_options(["--enable", "t1pks"])
        assert options.enabled_cluster_types == ["t1pks"]


class HmmDetectionTest(unittest.TestCase):
    def setUp(self):
        self.results_by_id = {
            "GENE_1": [
                FakeHSP("modelA", "GENE_1", 0, 10, 50, 0),
                FakeHSP("modelB", "GENE_1", 0, 10, 50, 0)
            ],
            "GENE_2": [
                FakeHSP("modelC", "GENE_2", 0, 10, 50, 0),
                FakeHSP("modelB", "GENE_2", 0, 10, 50, 0)
            ],
            "GENE_3": [
                FakeHSP("modelC", "GENE_3", 0, 10, 50, 0),
                FakeHSP("modelF", "GENE_3", 0, 10, 50, 0)
            ],
            "GENE_4": [
                FakeHSP("modelA", "GENE_4", 0, 10, 50, 0),
                FakeHSP("modelE", "GENE_4", 0, 10, 50, 0)
            ],
            "GENE_5": [
                FakeHSP("modelA", "GENE_5", 0, 10, 50, 0),
                FakeHSP("modelG", "GENE_5", 0, 10, 50, 0)
            ]
        }
        self.feature_by_id = {
            "GENE_1": DummyCDS(0, 30000, locus_tag="GENE_1"),
            "GENE_2": DummyCDS(30000, 50000, locus_tag="GENE_2"),
            "GENE_3": DummyCDS(70000, 90000, locus_tag="GENE_3"),
            "GENE_X": DummyCDS(95000, 100000, locus_tag="GENE_X"),  # no hits
            "GENE_4": DummyCDS(125000, 140000, locus_tag="GENE_4"),
            "GENE_5": DummyCDS(130000, 150000, locus_tag="GENE_5")
        }

        test_names = set(["modelA", "modelB", "modelC", "modelF", "modelG"])
        mock('signatures.get_signature_names', returns=test_names)

        self.rules = rule_parser.Parser([
                "MetaboliteA 10 5 modelA",
                "MetaboliteB 10 5 cds(modelA and modelB)",
                "MetaboliteC 10 5 (modelA and modelB)",
                "MetaboliteD 20 5 minimum(2,[modelC,modelB]) and modelA",
                "Metabolite0 1 3 modelF",
                "Metabolite1 1 3 modelG"]).rules
        self.features = []
        for gene_id in self.feature_by_id:
            self.features.append(self.feature_by_id[gene_id])
        self.features.sort(key=lambda x: x.location.start)  # vital for py3 < 3.5
        self.record = Record()
        self.record._record.seq = Seq("A"*150000)
        for feature in self.features:
            self.record.add_cds_feature(feature)

    def tearDown(self):
        # clear out any leftover config adjustments
        get_config().__dict__.clear()
        restore()

    def test_core(self):
        restore()  # don't mock signature names for this test

        def as_list_and_string(types, should_be_empty):
            " both list in options and ',' or ';' separated strings are fine "
            for case in [types, ",".join(types), ";".join(types)]:
                options = update_config({'enabled_cluster_types': case})
                if should_be_empty:
                    assert not core.check_options(options)
                else:
                    assert core.check_options(options)
                options.__dict__.clear()

        # should be no failing prerequisites
        assert core.check_prereqs() == []
        # always runs
        assert core.is_enabled(None)
        # default search types should be fine
        defaults = core.get_supported_cluster_types()
        as_list_and_string(defaults, True)
        # subsets of supported types should be fine
        as_list_and_string(defaults[:len(defaults)//2], True)
        # but adding an invalid cluster type should be caught
        as_list_and_string(defaults + ["bad_type"], False)

    def test_apply_cluster_rules(self):
        detected_types = hmm_detection.apply_cluster_rules(self.results_by_id,
                self.feature_by_id, self.rules)
        for gid in detected_types:
            detected_types[gid] = set(detected_types[gid].split("-"))
        expected_types = {
            "GENE_1": set(["MetaboliteA", "MetaboliteB", "MetaboliteC", "MetaboliteD"]),
            "GENE_2": set(["MetaboliteC", "MetaboliteD"]),
            "GENE_3": set(["Metabolite0"]),
            "GENE_4": set(["MetaboliteA"]),
            "GENE_5": set(["Metabolite1", "MetaboliteA"])
        }
        assert detected_types == expected_types

    def test_find_clusters(self):
        nseqdict = {"Metabolite0": "?", "Metabolite1": "?"}
        for i, gene_id in enumerate(self.feature_by_id):
            if gene_id == "GENE_X":
                continue
            clustertype = "Metabolite%d" % (i % 2)
            hmm_detection._update_sec_met_entry(self.feature_by_id[gene_id],
                             self.results_by_id[gene_id], clustertype, nseqdict)
        rules = {rule.name: rule for rule in self.rules}
        hmm_detection.find_clusters(self.record, rules)
        result_clusters = []
        for cluster in self.record.get_clusters():
            result_clusters.append(sorted(cds.get_name() for cds in cluster.cds_children))

        expected_clusters = [
            ["GENE_1", "GENE_2"],
            ["GENE_3"],
            ["GENE_4", "GENE_5"]
        ]
        assert result_clusters == expected_clusters

    def test_create_rules(self):
        restore()  # don't mock signature names for this test
        assert hmm_detection.create_rules([]) == []
        assert hmm_detection.create_rules(['nonexistant_cluster_type']) == []
        t1pks_rules = hmm_detection.create_rules(['t1pks'])
        assert len(t1pks_rules) == 1
        rule = t1pks_rules[0]
        assert rule.name == 't1pks'
        assert rule.cutoff == 20000
        assert rule.extent == 20000

    def test_get_overlaps_table(self):
        get_overlaps_table = hmm_detection.get_overlaps_table
        assert get_overlaps_table(DummyRecord()) == {}
        assert len(self.record.get_cds_features()) == 6
        assert get_overlaps_table(self.record) == {'GENE_1': 0, 'GENE_2': 1,
                                                   'GENE_3': 2, 'GENE_4': 4,
                                                   'GENE_5': 4, 'GENE_X': 3}

    def test_profiles_parsing(self):
        profiles = signatures.get_signature_profiles()
        assert len(profiles) == 220

    def test_store_detection(self):
        features = [DummyCluster(0, 10),
                    DummyFeature(10, 20),
                    DummyCluster(20, 30)]
        record = DummyRecord(features)
        assert len(record.get_clusters()) == 2
        for cluster in record.get_clusters():
            assert not cluster.detection_rules
        dummy_rule = rule_parser.DetectionRule("dummy", 10, 20, "a")
        hmm_detection.store_detection_details({"dummy": dummy_rule}, record)
        for feature in record.get_clusters():
            assert feature.products == ['dummy']
            assert feature.detection_rules == ['a']

    def test_store_detection_multitype(self):
        feature = DummyCluster(0, 50)
        feature.products = ["a", "b"]
        dummy_rule_a = rule_parser.DetectionRule("a", 10, 20, "c")
        dummy_rule_b = rule_parser.DetectionRule("b", 10, 20, "d")
        hmm_detection.store_detection_details({"a": dummy_rule_a, "b": dummy_rule_b},
                                              DummyRecord(features=[feature]))
        assert feature.detection_rules == ['c', 'd']
        assert feature.products == ['a', 'b']

    def test_filter(self):
        # fake HSPs all in one CDS with overlap > 20 and query_ids from the same equivalence group

        # not overlapping by > 20
        first = FakeHSP("AMP-binding", "A", 50, 90, 0.1, None)
        second = FakeHSP("A-OX", "A", 70, 100, 0.5, None)
        new, by_id = hmm_detection.filter_results([first, second], {"A": [first, second]})
        assert new == [first, second]
        assert by_id == {"A": [first, second]}

        # overlapping, in same group
        first.hit_end = 91
        assert hmm_detection.hsp_overlap_size(first, second) == 21
        new, by_id = hmm_detection.filter_results([first, second], {"A": [first, second]})
        assert new == [second]
        assert by_id == {"A": [second]}

        # overlapping, not in same group
        second.query_id = "none"
        new, by_id = hmm_detection.filter_results([first, second], {"A": [first, second]})
        assert new == [first, second]
        assert by_id == {"A": [first, second]}

        # not in the same CDS, but int he same group
        second.hit_id = "B"
        second.query_id = "A-OX"
        new, by_id = hmm_detection.filter_results([first, second], {"A": [first], "B": [second]})
        assert new == [first, second]
        assert by_id == {"A": [first], "B": [second]}

    def test_filter_multiple(self):
        # all in one CDS no overlap and the same query_ids -> cull all but the best score

        # not overlapping, not same query_id
        first = FakeHSP("AMP-binding", "A", 50, 60, 0.1, None)
        second = FakeHSP("A-OX", "A", 70, 100, 0.5, None)
        both = [first, second]
        by_id = {"A": [first, second]}
        new, by_id = hmm_detection.filter_result_multiple(list(both), dict(by_id))
        assert new == [first, second]
        assert by_id == {"A": [first, second]}

        # not overlapping, same query_id
        first.query_id = "A-OX"
        new, by_id = hmm_detection.filter_result_multiple(list(both), dict(by_id))
        assert new == [second]
        assert by_id == {"A": [second]}

        # not in same CDS, same query_id
        second.hit_id = "B"
        by_id = {"A": [first], "B": [second]}
        new, by_id = hmm_detection.filter_result_multiple(list(both), dict(by_id))
        assert new == [first, second]
        assert by_id == {"A": [first], "B": [second]}

    def test_filter_overlapping_genes(self):
        # overlapping CDS: take the longest CDS from within the same equivalence group
        def setup_and_run(results, features=None):
            if not features:
                features = []
                for hsp in results:
                    loc = FeatureLocation(hsp.hit_start, hsp.hit_end)
                    if hsp.query_id[0] == "A":
                        features.append(CDSFeature(loc, "dummy_trans", locus_tag=hsp.hit_id[0]))
                    else:
                        features.append(Feature(loc, feature_type="not_CDS"))
                        features[-1].locus_tag = hsp.hit_id[0]
            record = DummyRecord(features)
            overlaps = hmm_detection.get_overlaps_table(record)
            features_by_id = {i.locus_tag: i for i in features}
            results_by_id = {}
            for hsp in results:
                if hsp.hit_id not in results_by_id:
                    results_by_id[hsp.hit_id] = []
                results_by_id[hsp.hit_id].append(hsp)
            return hmm_detection.filter_result_overlapping_genes(
                    results, results_by_id, overlaps, features_by_id)

        # with overlap
        first = FakeHSP("AMP-binding", "A", 50, 71, 0.1, None)
        second = FakeHSP("A-OX", "B", 68, 100, 0.5, None)
        new_results, new_by_id = setup_and_run([first, second])
        assert new_results == [second]
        assert new_by_id == {"B": [second]}

        # without overlap
        first.hit_end = 67
        new_results, new_by_id = setup_and_run([first, second])
        assert new_results == [first, second]
        assert new_by_id == {"A": [first], "B": [second]}

        # make sure it catches a not-CDS if it snuck in
        with self.assertRaises(AssertionError):
            feature = DummyFeature(0, 20)
            feature.locus_tag = "none"
            setup_and_run([first], [feature])

    def test_equivalence_groups(self):
        group_file = path.get_full_path(os.path.dirname(__file__), "filterhmmdetails.txt")
        sets = []
        with open(group_file) as group_lines:
            sets = [set(line.strip().split(',')) for line in group_lines]

        # ensure they have at least two elements
        assert all(len(s) > 1 for s in sets)

        # ensure that the groups are disjoint
        for i, group in enumerate(sets):
            for other in sets[i + 1:]:
                assert group.isdisjoint(other)

    def test_hsp_overlap_size(self):
        overlap_size = hmm_detection.hsp_overlap_size
        first = FakeHSP("A", "A", 50, 60, 0., None)
        second = FakeHSP("B", "B", 70, 100, 0., None)
        # no overlap
        assert overlap_size(first, second) == 0
        first.hit_end = 70
        # still no overlap, end isn't inclusive
        assert overlap_size(first, second) == 0
        # a mix of second starting inside first
        for i in range(1, 30):
            first.hit_end += 1
            assert overlap_size(first, second) == i
        # second wholly contained
        first.hit_end = 110
        assert overlap_size(first, second) == 30

        # first inside second
        first.hit_start = 75
        assert overlap_size(first, second) == 25

        # first inside second, but direction reversed
        first.hit_end = 50
        with self.assertRaises(AssertionError):
            overlap_size(first, second)
