# License: GNU Affero General Public License v3 or later
# A copy of GNU AGPL v3 should have been included in this software package in LICENSE.txt.

# for test files, silence irrelevant and noisy pylint warnings
# pylint: disable=no-self-use,protected-access,missing-docstring

import glob
import os
import unittest

from Bio.Seq import Seq

from antismash.common import path
from antismash.common.hmm_rule_parser import rule_parser, cluster_prediction as hmm_detection  # TODO: redo tests
from antismash.common.hmm_rule_parser.test.helpers import check_hmm_signatures
from antismash.common.secmet import Record
from antismash.common.test.helpers import DummyCDS, FakeHSPHit
from antismash.config import get_config
import antismash.detection.hmm_detection as core
from antismash.detection.hmm_detection import signatures


class HmmDetectionTest(unittest.TestCase):
    def setUp(self):
        self.rules_file = path.get_full_path(__file__, "..", "cluster_rules.txt")
        self.signature_file = path.get_full_path(__file__, "..", "data", "hmmdetails.txt")
        self.signature_names = {sig.name for sig in core.get_signature_profiles()}
        self.filter_file = path.get_full_path(__file__, "..", "filterhmmdetails.txt")
        self.results_by_id = {
            "GENE_1": [
                FakeHSPHit("modelA", "GENE_1", 0, 10, 50, 0),
                FakeHSPHit("modelB", "GENE_1", 0, 10, 50, 0)
            ],
            "GENE_2": [
                FakeHSPHit("modelC", "GENE_2", 0, 10, 50, 0),
                FakeHSPHit("modelB", "GENE_2", 0, 10, 50, 0)
            ],
            "GENE_3": [
                FakeHSPHit("modelC", "GENE_3", 0, 10, 50, 0),
                FakeHSPHit("modelF", "GENE_3", 0, 10, 50, 0)
            ],
            "GENE_4": [
                FakeHSPHit("modelA", "GENE_4", 0, 10, 50, 0),
                FakeHSPHit("modelE", "GENE_4", 0, 10, 50, 0)
            ],
            "GENE_5": [
                FakeHSPHit("modelA", "GENE_5", 0, 10, 50, 0),
                FakeHSPHit("modelG", "GENE_5", 0, 10, 50, 0)
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

        self.test_names = {"modelA", "modelB", "modelC", "modelF", "modelG",
                           "a", "b", "c", "d"}

        self.rules = rule_parser.Parser("\n".join([
                "RULE MetaboliteA CUTOFF 10 EXTENT 5 CONDITIONS modelA",
                "RULE MetaboliteB CUTOFF 10 EXTENT 5 CONDITIONS cds(modelA and modelB)",
                "RULE MetaboliteC CUTOFF 10 EXTENT 5 CONDITIONS (modelA and modelB)",
                "RULE MetaboliteD CUTOFF 20 EXTENT 5 CONDITIONS minimum(2,[modelC,modelB]) and modelA",
                "RULE Metabolite0 CUTOFF 1 EXTENT 3 CONDITIONS modelF",
                "RULE Metabolite1 CUTOFF 1 EXTENT 3 CONDITIONS modelG"]), self.test_names).rules
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

    def test_overlaps_but_not_contains(self):
        # should get gene2 and gene3
        rules = rule_parser.Parser("\n".join([
                "RULE Overlap CUTOFF 25 EXTENT 5 CONDITIONS modelB and modelF "
                "RULE OverlapImpossible CUTOFF 25 EXTENT 5 CONDITIONS modelA and modelF"]), self.test_names).rules
        detected_types, cluster_type_hits = hmm_detection.apply_cluster_rules(self.record, self.results_by_id, rules)
        assert detected_types == {"GENE_2": {"Overlap": {"modelB"}},
                                  "GENE_3": {"Overlap": {"modelF"}}}

        assert cluster_type_hits == {"Overlap": {"GENE_2", "GENE_3"}}

        # only 1 cluster should be found, since it requires both genes
        # if forming clusters by .is_contained_by(), 2 clusters will be formed
        # if finding rule hits uses .is_contained_by(), no clusters will be formed
        rules_by_name = {rule.name: rule for rule in rules}
        clusters = hmm_detection.find_clusters(self.record, cluster_type_hits, rules_by_name)
        assert len(clusters) == 1
        assert clusters[0].product == "Overlap"
        assert clusters[0].core_location.start == 30000
        assert clusters[0].core_location.end == 90000

    def test_core(self):
        # should be no failing prerequisites
        assert core.check_prereqs() == []
        # always runs
        assert core.is_enabled(None)

    def test_apply_cluster_rules(self):
        detected_types, cluster_type_hits = hmm_detection.apply_cluster_rules(self.record, self.results_by_id,
                                                                              self.rules)
        for gid in detected_types:
            detected_types[gid] = set(detected_types[gid])
        expected_types = {
            "GENE_1": set(["MetaboliteA", "MetaboliteB", "MetaboliteC", "MetaboliteD"]),
            "GENE_2": set(["MetaboliteC", "MetaboliteD"]),
            "GENE_3": set(["Metabolite0"]),
            "GENE_4": set(["MetaboliteA"]),
            "GENE_5": set(["Metabolite1", "MetaboliteA"])
        }
        assert detected_types == expected_types

        assert cluster_type_hits == {"MetaboliteA": {"GENE_1", "GENE_4", "GENE_5"},
                                     "MetaboliteB": {"GENE_1"},
                                     "MetaboliteC": {"GENE_1", "GENE_2"},
                                     'MetaboliteD': {'GENE_1', 'GENE_2'},
                                     'Metabolite0': {'GENE_3'},
                                     'Metabolite1': {'GENE_5'}}

    def test_find_clusters(self):
        cds_features_by_type = {"MetaboliteA": {"GENE_1", "GENE_4", "GENE_5"},
                                "MetaboliteB": {"GENE_1"},
                                "MetaboliteC": {"GENE_1", "GENE_2"},
                                'MetaboliteD': {'GENE_1', 'GENE_2'},
                                'Metabolite0': {'GENE_3'},
                                'Metabolite1': {'GENE_5'}}
        rules = {rule.name: rule for rule in self.rules}
        for cluster in hmm_detection.find_clusters(self.record, cds_features_by_type, rules):
            self.record.add_cluster(cluster)
        assert len(self.record.get_clusters()) == 7
        cluster_products = sorted([cluster.product for cluster in self.record.get_clusters()])
        assert cluster_products == sorted(["Metabolite%s" % i for i in "01AABCD"])
        self.record.create_superclusters()
        assert len(self.record.get_superclusters()) == 3
        self.record.create_regions()
        assert len(self.record.get_regions()) == 3
        result_regions = []
        for region in self.record.get_regions():
            result_regions.append(sorted(cds.get_name() for cds in region.cds_children))

        expected_regions = [
            ["GENE_1", "GENE_2"],
            ["GENE_3"],
            ["GENE_4", "GENE_5"]
        ]
        assert result_regions == expected_regions

    def test_create_rules(self):
        rules = hmm_detection.create_rules(self.rules_file, self.signature_names)
        assert len(rules) == 47
        t1pks_rules = [rule for rule in rules if rule.name == "t1pks"]
        assert len(t1pks_rules) == 1
        rule = t1pks_rules[0]
        assert rule.name == 't1pks'
        assert rule.cutoff == 20000
        assert rule.extent == 20000

    def test_profiles_parsing(self):
        profiles = signatures.get_signature_profiles()
        assert len(profiles) == 226  # ensures we don't delete any by accident

    def test_filter(self):
        # fake HSPs all in one CDS with overlap > 20 and query_ids from the same equivalence group

        # not overlapping by > 20
        first = FakeHSPHit("AMP-binding", "A", 50, 90, 0.1, None)
        second = FakeHSPHit("A-OX", "A", 70, 100, 0.5, None)
        new, by_id = hmm_detection.filter_results([first, second], {"A": [first, second]},
                                                  self.filter_file, self.signature_names)
        assert new == [first, second]
        assert by_id == {"A": [first, second]}

        # overlapping, in same group
        first.hit_end = 91
        assert hmm_detection.hsp_overlap_size(first, second) == 21
        new, by_id = hmm_detection.filter_results([first, second], {"A": [first, second]},
                                                  self.filter_file, self.signature_names)
        assert new == [second]
        assert by_id == {"A": [second]}

        # overlapping, not in same group
        second.query_id = "none"
        new, by_id = hmm_detection.filter_results([first, second], {"A": [first, second]},
                                                  self.filter_file, self.signature_names)
        assert new == [first, second]
        assert by_id == {"A": [first, second]}

        # not in the same CDS, but int he same group
        second.hit_id = "B"
        second.query_id = "A-OX"
        new, by_id = hmm_detection.filter_results([first, second], {"A": [first], "B": [second]},
                                                  self.filter_file, self.signature_names)
        assert new == [first, second]
        assert by_id == {"A": [first], "B": [second]}

    def test_filter_multiple(self):
        # all in one CDS no overlap and the same query_ids -> cull all but the best score

        # not overlapping, not same query_id
        first = FakeHSPHit("AMP-binding", "A", 50, 60, 0.1, None)
        second = FakeHSPHit("A-OX", "A", 70, 100, 0.5, None)
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
        first = FakeHSPHit("A", "A", 50, 60, 0., None)
        second = FakeHSPHit("B", "B", 70, 100, 0., None)
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

    def test_hmm_files_and_details_match(self):
        data_dir = path.get_full_path(os.path.dirname(__file__), "data", "")
        details_files = {prof.path for prof in signatures.get_signature_profiles()}
        details_files = {filepath.replace(data_dir, "") for filepath in details_files}
        data_dir_contents = set(glob.glob(data_dir + "*.hmm"))
        data_dir_contents = {filepath.replace(data_dir, "") for filepath in data_dir_contents}
        # ignore bgc_seeds.hmm for the sake of comparison, it's a generated aggregate
        data_dir_contents.discard("bgc_seeds.hmm")
        missing_files = details_files - data_dir_contents
        assert not missing_files
        extra_files = data_dir_contents - details_files
        assert not extra_files
        # finally, just to be sure
        assert data_dir_contents == details_files


class TestSignatureFile(unittest.TestCase):
    def test_details(self):
        data_dir = path.get_full_path(os.path.dirname(__file__), 'data')
        check_hmm_signatures(os.path.join(data_dir, 'hmmdetails.txt'), data_dir)
