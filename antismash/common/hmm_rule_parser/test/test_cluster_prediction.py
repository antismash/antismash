# License: GNU Affero General Public License v3 or later
# A copy of GNU AGPL v3 should have been included in this software package in LICENSE.txt.

# for test files, silence irrelevant and noisy pylint warnings
# pylint: disable=use-implicit-booleaness-not-comparison,protected-access,missing-docstring,too-many-public-methods

from dataclasses import FrozenInstanceError
import unittest
from unittest.mock import patch

from antismash.common.hmm_rule_parser import cluster_prediction, rule_parser, structures
from antismash.common.secmet.features import Protocluster
from antismash.common.secmet.locations import CompoundLocation, FeatureLocation
from antismash.common.secmet.qualifiers.gene_functions import GeneFunction
from antismash.common.test.helpers import DummyRecord, DummyCDS, DummyProtocluster, FakeHSPHit

from .helpers import create_ruleset


class DummyConditions(rule_parser.Conditions):
    """ so a DetectionRule can be created without failing its internal checks """
    def __init__(self):
        super().__init__(negated=False)

    def contains_positive_condition(self):
        return True

# NOTE: the rest of the cluster_prediction tests are still in hmm_detection tests


class TestAncillary(unittest.TestCase):
    def test_ancillary_on_cutoff_boundary(self):
        features_by_id = {
            "a": DummyCDS(3082, 6555, locus_tag="a"),
            "b": DummyCDS(7650, 9989, locus_tag="b"),
            "c": DummyCDS(10000, 12618, locus_tag="c"),
            "d": DummyCDS(13840, 15411, locus_tag="d"),
            "e": DummyCDS(15596, 17275, locus_tag="e"),
        }
        features = list(features_by_id.values())
        record = DummyRecord(features)

        results_by_id = {
            name: [FakeHSPHit(name.upper(), name, 0, 10, 50, 0)] for name in features_by_id
        }
        signature_names = set(name.upper() for name in features_by_id)
        rule_text = "\n".join([
            "RULE rule_name",
            "CATEGORY category",
            "CUTOFF 5",
            "NEIGHBOURHOOD 10",
            "CONDITIONS A and B and C and D and E",
        ])
        rules = rule_parser.Parser(rule_text, signature_names, {"category"}).rules

        cds_to_rules, rules_to_cds = cluster_prediction.apply_cluster_rules(record, results_by_id, rules)

        # with a cutoff of 5, only 'c' is in range of everything else, and if
        # the other genes that help satisfy the conditions aren't included in 
        # the ancillary hits, then the resulting protocluster will not be long
        # enough to cover all the genes that were required to satisfy the conditions
        assert len(rules_to_cds) == 1  # only one rule exists
        # all genes should be core genes for the rule
        assert rules_to_cds[rules[0].name] == set(features_by_id)
        # inversely, all genes should have links back to that rule
        assert sorted(cds_to_rules) == sorted(features_by_id)


class TestRedundancy(unittest.TestCase):
    def setUp(self):
        self.record = DummyRecord()
        for cds in [DummyCDS(start=50, end=80), DummyCDS(start=110, end=140), DummyCDS(start=150, end=180)]:
            self.record.add_cds_feature(cds)
        superior = rule_parser.DetectionRule("superior", "category", 10, 10, DummyConditions())
        inferior = rule_parser.DetectionRule("inferior", "category", 10, 10, DummyConditions(), superiors=["superior"])
        irrelevant = rule_parser.DetectionRule("irrelevant", "category", 10, 10, DummyConditions())
        self.rules_by_name = {rule.name: rule for rule in [superior, inferior, irrelevant]}

    def remove(self, clusters):
        return cluster_prediction.remove_redundant_protoclusters(clusters, self.rules_by_name, self.record)

    def create_cluster(self, rule_name, start, end):
        rule = self.rules_by_name[rule_name]
        core = FeatureLocation(start, end)
        surrounds = FeatureLocation(max(0, start - rule.neighbourhood), end + rule.neighbourhood)
        return Protocluster(core, surrounds, tool="testing", cutoff=rule.cutoff,
                            neighbourhood_range=rule.neighbourhood, product=rule_name,
                            detection_rule="rule text")

    def test_alone(self):
        clusters = [self.create_cluster("inferior", 50, 140)]
        assert clusters == self.remove(clusters)

    def test_non_overlap(self):
        clusters = [self.create_cluster("inferior", 50, 140),
                    self.create_cluster("superior", 150, 180)]
        assert clusters == self.remove(clusters)

    def test_not_relevant_equal(self):
        clusters = [self.create_cluster("inferior", 50, 140),
                    self.create_cluster("irrelevant", 50, 140)]
        assert clusters == self.remove(clusters)

    def test_not_relevant_contained(self):
        clusters = [self.create_cluster("inferior", 110, 140),
                    self.create_cluster("irrelevant", 50, 180)]
        assert clusters == self.remove(clusters)

    def test_not_relevant_larger(self):
        clusters = [self.create_cluster("inferior", 50, 180),
                    self.create_cluster("irrelevant", 110, 140)]
        assert clusters == self.remove(clusters)

    def test_contained(self):
        clusters = [self.create_cluster("inferior", 110, 140),
                    self.create_cluster("superior", 50, 180)]
        assert self.remove(clusters) == [clusters[1]]

    def test_equal(self):
        clusters = [self.create_cluster("inferior", 110, 140),
                    self.create_cluster("superior", 110, 140)]
        assert self.remove(clusters) == [clusters[1]]

    def test_larger(self):
        clusters = [self.create_cluster("inferior", 50, 180),
                    self.create_cluster("superior", 110, 140)]
        assert self.remove(clusters) == [clusters[1]]

    def test_adjacent_with_overlap(self):
        # no intersection of core genes at all
        # one gene just overlaps slightly with the previous
        existing_end = self.record.get_cds_features()[-1].location.end
        self.record.add_cds_feature(DummyCDS(start=existing_end - 10, end=existing_end + 20))
        clusters = [
            self.create_cluster("superior", 0, existing_end),
            self.create_cluster("inferior", existing_end - 10, existing_end + 20),
        ]
        # that adjacent cluster should not be discarded as redundant
        assert self.remove(clusters) == clusters

    def test_neighbourhoods_dont_matter(self):
        neighbourhood = self.rules_by_name["superior"].neighbourhood
        for new_neighbourhood in [neighbourhood - 10, neighbourhood + 10]:
            self.rules_by_name["inferior"].neighbourhood = new_neighbourhood
            self.test_larger()
            self.test_equal()
            self.test_contained()

    def test_cutoffs_dont_matter(self):
        cutoff = self.rules_by_name["superior"].cutoff
        for new_cutoff in [cutoff - 10, cutoff + 10]:
            self.rules_by_name["inferior"].cutoff = new_cutoff
            self.test_larger()
            self.test_equal()
            self.test_contained()


class TestDynamic(unittest.TestCase):
    def test_find_dynamic(self):
        expected_a = {"cds_name": [structures.DynamicHit("prof_a", "cds_name")]}
        expected_b = {"cds_name": [structures.DynamicHit("prof_b", "cds_name")]}
        profile_a = structures.DynamicProfile("prof_a", "desc a", lambda record, _: expected_a)
        profile_b = structures.DynamicProfile("prof_b", "desc b", lambda record, _: expected_b)
        results = cluster_prediction.find_dynamic_hits(DummyRecord(), [profile_a, profile_b], {})
        assert results["cds_name"] == expected_a["cds_name"] + expected_b["cds_name"]

    @patch.object(cluster_prediction, "find_hmmer_hits", return_value={})
    @patch.object(cluster_prediction, "get_signature_profiles", return_value={})
    def test_full(self, _patched_find, _patched_sigs):
        cdses = [DummyCDS(locus_tag="A", start=0, end=6), DummyCDS(locus_tag="B", start=6, end=12)]
        record = DummyRecord(features=cdses)

        # create a dummy dynamic profile
        def find_a(rec, _hmmer_hits):
            hits = {}
            for cds in rec.get_cds_features():
                if cds.get_name() == "A":
                    hits[cds.get_name()] = [structures.DynamicHit(cds.get_name(), "a_finder")]
            return hits
        profile = structures.DynamicProfile("a_finder", "desc", find_a)
        # make sure the 'profile' functions as expected
        assert cdses[0].get_name() in profile.find_hits(record, {})

        # build a dummy rule that will search for this hit
        condition = rule_parser.SingleCondition(False, "a_finder")
        rule = rule_parser.DetectionRule("test-name", "Other", 5000, 5000, condition)
        ruleset = create_ruleset([rule], categories=set("Other"), dynamic_profiles={profile.name: profile})
        results = cluster_prediction.detect_protoclusters_and_signatures(record, ruleset, "test_tool")
        assert results
        assert results.cds_by_cluster
        assert results.protoclusters
        proto = results.protoclusters[0]
        assert proto.product == "test-name"

        results.annotate_cds_features()
        assert cdses[0].sec_met.domains[0].name == "a_finder"


class TestRecordContainment(unittest.TestCase):
    def test_neighbourhood(self):
        simple_hits = {"cds_name": [structures.DynamicHit("cds_name", "profile")]}
        profile = structures.DynamicProfile("profile", "desc", lambda record, _: simple_hits)

        cds = DummyCDS(locus_tag="cds_name", start=10, end=24)
        record = DummyRecord(seq="A" * 30, features=[cds])

        rule_text = "RULE A CATEGORY Cat CUTOFF 10 NEIGHBOURHOOD 5 CONDITIONS profile"
        rules = rule_parser.Parser(rule_text, {"profile"}, {"Cat"}).rules
        rules[0].neighbourhood = 8  # not KB, just bases, and enough to extend over the end
        ruleset = create_ruleset(tuple(rules), dynamic_profiles={"profile": profile})

        results = cluster_prediction.detect_protoclusters_and_signatures(record, ruleset)
        assert results.protoclusters
        for protocluster in results.protoclusters:
            assert protocluster.location.start == 2
            assert protocluster.location.end == len(record)
            assert protocluster.core_location == cds.location

        record.make_circular()

        results = cluster_prediction.detect_protoclusters_and_signatures(record, ruleset)
        assert results.protoclusters
        for protocluster in results.protoclusters:
            assert len(protocluster.location.parts) == 2
            assert protocluster.location.parts[0].start == 2
            assert protocluster.location.parts[0].end == len(record)
            assert protocluster.location.parts[1].start == 0
            assert protocluster.location.parts[1].end == 2
            assert protocluster.core_location == cds.location


class TestMultipliers(unittest.TestCase):
    def test_create_rules(self):
        text = "RULE A CATEGORY Cat CUTOFF 10 NEIGHBOURHOOD 5 CONDITIONS A"
        # with default multipliers
        with patch("builtins.open", unittest.mock.mock_open(read_data=text)):
            rule = cluster_prediction.create_rules(["dummy.file"], {"A"}, {"Cat"}, {})[0]
        assert rule.cutoff == 10_000
        assert rule.neighbourhood == 5_000

        # with custom multipliers
        multipliers = cluster_prediction.Multipliers(
            cutoff=1.5,
            neighbourhood=0.5,
        )
        with patch("builtins.open", unittest.mock.mock_open(read_data=text)):
            multiplied = cluster_prediction.create_rules(["dummy.file"], {"A"}, {"Cat"},
                                                         multipliers=multipliers,
                                                         )[0]
        # make sure the multipliers were used
        assert multiplied.cutoff == rule.cutoff * multipliers.cutoff
        assert multiplied.neighbourhood == rule.neighbourhood * multipliers.neighbourhood

    def test_multipliers_unchangeable(self):
        multipliers = structures.Multipliers(1, 5)
        assert multipliers.cutoff == 1
        assert multipliers.neighbourhood == 5
        with self.assertRaises(FrozenInstanceError):
            multipliers.cutoff = 1
        with self.assertRaises(FrozenInstanceError):
            multipliers.neighbourhood = 5


class TestDomainAnnotations(unittest.TestCase):
    def test_inferiors_not_annotated(self):
        cds_a = DummyCDS(start=10, end=40, locus_tag="a")
        cds_b = DummyCDS(start=50, end=80, locus_tag="b")
        record = DummyRecord(features=[cds_a, cds_b])
        clusters = [DummyProtocluster(start=10, end=80, product="superior")]  # inferior cluster already discarded
        hsps = {
            "a": [FakeHSPHit("sup1", "a"), FakeHSPHit("inf1", "a")],
            "b": [FakeHSPHit("inf2", "a")],
        }
        domains_by_cluster = {
            "a": {
                "superior": {"sup1"},
                "inferior": {"inf1"},
            },
            "b": {
                "inferior": {"inf2"},
            },
        }

        results = cluster_prediction.build_results(clusters, record, "dummy_tool", hsps, domains_by_cluster,
                                                   True, structures.Multipliers())
        # neither gene should have function annotations yet
        assert len(cds_a.gene_functions) == 0
        assert len(cds_b.gene_functions) == 0

        results.annotate_cds_features()
        # now A should have the superior rule's functions as core and inferiors as additional
        assert len(cds_a.gene_functions) == 2
        assert cds_a.gene_functions.get_by_function(GeneFunction.CORE)
        assert cds_a.gene_functions.get_by_function(GeneFunction.ADDITIONAL)
        # and B should have the inferiors as additional
        assert len(cds_b.gene_functions) == 1
        assert not cds_b.gene_functions.get_by_function(GeneFunction.CORE)
        assert cds_b.gene_functions.get_by_function(GeneFunction.ADDITIONAL)


class TestCircularity(unittest.TestCase):
    def setUp(self):
        text = (
            "RULE A CATEGORY Cat CUTOFF 1 NEIGHBOURHOOD 1 CONDITIONS profA and profB "
            "RULE B CATEGORY Cat CUTOFF 0 NEIGHBOURHOOD 0 CONDITIONS profC "
        )
        rules = rule_parser.Parser(text, {"profA", "profB", "profC"}, {"Cat"}).rules
        self.rules = {rule.name: rule for rule in rules}
        self.rule = rules[0]
        self.rules["A"].neighbourhood = 100  # the sequence can then be much smaller
        self.rules["A"].cutoff = 100
        self.hits = {"A": {"left_core", "right_core"}, "B": {"left_distant", "right_distant"}}
        self.find_func = cluster_prediction.find_protoclusters
        cdses = [
            DummyCDS(90, 93, 1, locus_tag="left_distant"),
            DummyCDS(195, 198, 1, locus_tag="left_neighbour"),
            DummyCDS(200, 206, 1, locus_tag="left_core"),
            DummyCDS(212, 215, 1, locus_tag="right_core"),
            DummyCDS(220, 223, 1, locus_tag="right_neighbour"),
            DummyCDS(350, 356, 1, locus_tag="right_distant"),
        ]
        self.length = 800
        self.cdses = {cds.locus_tag: cds for cds in cdses}
        self.record = DummyRecord(features=list(self.cdses.values()), seq="A" * self.length, circular=True)

    def check(self):
        assert self.record.is_circular()
        results = self.find_func(self.record, self.hits, self.rules, {}, {})
        assert len(results) == 3

        proto = [p for p in results if p.product == "A"][0]
        print(proto, proto.core_location)
        # the "distant" CDSes shouldn't be included
        included = sorted(self.record.get_cds_features_within_location(proto.location))
        expected = sorted(cds for cds in self.cdses.values() if "distant" not in cds.get_name())
        assert included == expected

        # the core location should match the core genes only
        cores = [cds.get_name() for cds in self.record.get_cds_features_within_location(proto.core_location,
                                                                                        with_overlapping=False)]
        assert cores == ["left_core", "right_core"]
        assert 1 <= len(proto.core_location.parts) <= 2
        assert proto.core_location.parts[0].start == self.cdses["left_core"].location.parts[0].start
        assert proto.core_location.parts[-1].end == self.cdses["right_core"].location.parts[-1].end

        # the neighbourhood should match the extension of the core
        assert str(proto.location) == str(self.record.extend_location(proto.core_location, self.rule.neighbourhood))

        # the neighbourhood must have multiple parts if the core does, but may have an extra
        assert len(proto.location.parts) >= len(proto.core_location.parts)

        return proto

    def test_merging_uses_first_and_last(self):
        # unlike linear merges, the circular needs to check the first and last
        # to see if they're close enough over the origin, but if the first and
        # last aren't correctly selected for distance checks, then the merge
        # will do very bad things
        cdses = []
        locations = [FeatureLocation(start, end, 1) for start, end in [
            (3091, 3092), (3105, 3111), (6882, 6885), (6887, 6889)
        ]]
        cdses = [DummyCDS(location.start, location.end) for location in locations]
        record = DummyRecord(seq="A", features=cdses, circular=True, length=max(loc.end for loc in locations))
        hits = {"A": set(cds.get_name() for cds in cdses)}
        self.rules["A"].cutoff = 20
        self.rules["A"].neighbourhood = 0
        results = self.find_func(record, hits, self.rules, {}, {})
        # both first two and last two should have been merged
        assert len(results) == len(cdses) - 2
        assert results[0].location == FeatureLocation(locations[0].start, locations[1].end, 1)
        assert results[1].location == FeatureLocation(locations[2].start, locations[3].end, 1)

    def test_without_crossing_origin(self):
        proto = self.check()
        assert len(proto.location.parts) == 1

    def test_crossing_origin_between_cores(self):
        self.record.rotate(209)
        assert self.cdses["left_core"] > self.cdses["right_core"]
        proto = self.check()
        assert len(proto.core_location.parts) == 2
        assert len(proto.location.parts) == 2

    def test_left_neighbourhood_beyond_origin(self):
        self.record.rotate(199)
        assert self.cdses["left_neighbour"] > self.cdses["left_core"]
        proto = self.check()
        assert len(proto.core_location.parts) == 1
        assert len(proto.location.parts) == 2

    def test_right_neighbourhood_beyond_origin(self):
        self.record.rotate(216)
        assert self.cdses["left_neighbour"] > self.cdses["right_neighbour"]
        proto = self.check()
        assert len(proto.core_location.parts) == 1
        assert len(proto.location.parts) == 2

    def test_splitting_core_gene(self):
        for target in ["left_core", "right_core"]:
            location = self.cdses[target].location
            self.record.rotate((location.start + location.end) // 2)
            assert len(self.cdses[target].location.parts) > 1
            proto = self.check()
            assert len(proto.core_location.parts) == 2
            assert len(proto.location.parts) == 2

    def test_overlapping_cores_at_origin(self):
        self.record = DummyRecord(length=7663439, circular=True)
        cores = [
            DummyCDS(start=7633885, end=7644192, strand=1, locus_tag="pre_origin"),
            DummyCDS(location=CompoundLocation([
                FeatureLocation(7644229, len(self.record), 1),
                FeatureLocation(0, 19, 1),
            ]), locus_tag="mostly_pre_origin"),
            DummyCDS(location=CompoundLocation([
                FeatureLocation(7663427, len(self.record), 1),
                FeatureLocation(0, 9008, 1),
            ]), locus_tag="mostly_post_origin"),
        ]
        self.cdses = {cds.get_name(): cds for cds in cores}
        self.cdses["unrelated"] = DummyCDS(location=FeatureLocation(2100, 21000, 1), locus_tag="unrelated")
        for cds in self.cdses.values():
            self.record.add_cds_feature(cds)
        self.hits = {"A": set(core.get_name() for core in cores)}
        self.rules["A"].cutoff = 100
        protos = self.find_func(self.record, self.hits, self.rules, {}, {})
        assert len(protos) == 1
        proto = protos[0]
        assert proto.core_location == CompoundLocation([
            FeatureLocation(cores[0].start, len(self.record), 1),
            FeatureLocation(0, cores[-1].end, 1),
        ])
        for cds in self.cdses.values():
            if cds.get_name() == "unrelated":
                assert not cds.is_contained_by(proto.core_location)
            else:
                assert cds.is_contained_by(proto.core_location)


class TestLocationExtension(unittest.TestCase):
    def setUp(self):
        self.func = cluster_prediction._extend_area_location
        self.record = DummyRecord(seq="A" * 100)

    def test_overlapping_wraparound(self):
        expected = FeatureLocation(0, len(self.record), 1)
        result = self.func(FeatureLocation(28, 70, 1), 40, self.record)
        assert result == expected, f"{str(result)=} != {str(expected)=}"

    def test_strands_forward(self):
        expected = FeatureLocation(40, 70, 1)
        result = self.func(FeatureLocation(50, 60, -1), 10, self.record)
        assert expected == result

    def test_trivial(self):
        assert not self.record.is_circular()
        location = FeatureLocation(50, 60, 1)
        distance = 30
        result = self.func(location, distance, self.record)
        assert len(result.parts) == len(location.parts)
        assert result.start == location.start - distance
        assert result.end == location.end + distance
        assert len(result) == len(location) + distance * 2

    def test_non_circular_compound(self):
        assert not self.record.is_circular()
        parts = [
            FeatureLocation(20, 30, 1),
            FeatureLocation(50, 60, 1),
        ]
        location = CompoundLocation(parts)
        with self.assertRaisesRegex(ValueError, "too many sub-locations"):
            self.func(location, 10, self.record)

    def test_forward_crossing(self):
        self.record.make_circular()
        parts = [
            FeatureLocation(50, 100, 1),
            FeatureLocation(0, 30, 1),
        ]
        location = CompoundLocation(parts)
        distance = 5
        result = self.func(location, distance, self.record)
        assert len(result) == len(location) + distance * 2
        assert result.parts[0].start == parts[0].start - distance
        assert result.parts[0].end == parts[0].end == len(self.record)
        assert result.parts[1].start == parts[1].start == 0
        assert result.parts[1].end == parts[1].end + distance
