# License: GNU Affero General Public License v3 or later
# A copy of GNU AGPL v3 should have been included in this software package in LICENSE.txt.

# for test files, silence irrelevant and noisy pylint warnings
# pylint: disable=use-implicit-booleaness-not-comparison,protected-access,missing-docstring

""" Tests for building predictions from the rule detection results

    The rule detection results hold, for each protocluster, a CDSResults per
    CDS with hits inside the protocluster. Each of those has:
     - domains: every domain hit found in the CDS
     - definition_domains: for each rule, the domains of the CDS that fired
       that rule, which can be empty

    The same CDS can lie in multiple protoclusters, so it can have an entry for
    each of those rules. In the tests, THIS_RULE is the rule a prediction is
    built for and OTHER_RULE is a rule of an overlapping subcluster.
"""

from types import SimpleNamespace
import unittest
from unittest.mock import patch, sentinel

from antismash.common.hmm_rule_parser.cluster_prediction import CDSResults
from antismash.common.secmet.locations import CompoundLocation, FeatureLocation
from antismash.common.secmet.qualifiers import SecMetQualifier
from antismash.common.secmet.test.helpers import DummyCDS, DummyProtocluster
from antismash.detection.subclusters import predictions
from antismash.detection.subclusters.predictions import (
    CDSDomainHit,
    SubclusterPrediction,
    _build_domain_hits,
    build_predictions,
)
from antismash.detection.subclusters.signatures import SubclusterHmmSignature

THIS_RULE = "this_rule"
OTHER_RULE = "other_rule"

# signatures A to D, of which D isn't used by any rule
SIGNATURES = {
    name: SubclusterHmmSignature(name, f"{name} description", 10, f"{name}.hmm",
                                 accession=f"PF0000{i}.1")
    for i, name in enumerate(["A", "B", "C", "D"], start=1)
}


def build_domain(name, bitscore=50, evalue=1e-10):
    return SecMetQualifier.Domain(name, evalue, bitscore, 1, "test")


def build_cds_result(cds_name, hits, fired):
    """ Builds the results for a single CDS

        Arguments:
            cds_name: the name of the CDS
            hits: the names of the domains hit in the CDS, or the domains themselves
            fired: a dictionary mapping rule name to the domains that fired that rule
    """
    domains = [build_domain(hit) if isinstance(hit, str) else hit for hit in hits]
    return CDSResults(DummyCDS(locus_tag=cds_name), domains, fired)


def build_rule(name, conditions="(A and B)"):
    return SimpleNamespace(name=name, conditions=conditions)


def cds_and_domain_names(hits):
    """ Reduces domain hits to (CDS name, domain name) pairs, for simpler comparisons """
    return [(hit.cds_name, hit.domain_name) for hit in hits]


class TestBuildDomainHits(unittest.TestCase):
    def setUp(self):
        signature_patch = patch.object(predictions, "get_signatures", return_value=SIGNATURES)
        signature_patch.start()
        self.addCleanup(signature_patch.stop)

    def test_hit_details(self):
        domain = build_domain("A", bitscore=50, evalue=1e-20)
        cds_result = build_cds_result("cds1", hits=[domain], fired={THIS_RULE: {"A"}})

        # the details of the signature are combined with the scores of the hit
        assert _build_domain_hits(THIS_RULE, [cds_result]) == [
            CDSDomainHit(domain_name="A", domain_description="A description",
                         domain_accession="PF00001.1", cds_name="cds1",
                         evalue=1e-20, bitscore=50),
        ]

    def test_only_fired_domains(self):
        cds_results = [
            # C only fired for an overlapping subcluster
            build_cds_result("cds1", hits=["A", "B", "C"],
                             fired={THIS_RULE: {"A", "B"}, OTHER_RULE: {"C"}}),
            # only part of the overlapping subcluster
            build_cds_result("cds2", hits=["C"], fired={THIS_RULE: set(), OTHER_RULE: {"C"}}),
            # D isn't used by any rule
            build_cds_result("cds3", hits=["D"], fired={THIS_RULE: set()}),
        ]
        hits = _build_domain_hits(THIS_RULE, cds_results)
        assert cds_and_domain_names(hits) == [("cds1", "A"), ("cds1", "B")]

    def test_multiple_cdses(self):
        cds_results = [
            build_cds_result("cds1", hits=["A"], fired={THIS_RULE: {"A"}}),
            build_cds_result("cds2", hits=["A", "B"], fired={THIS_RULE: {"A", "B"}}),
        ]
        hits = _build_domain_hits(THIS_RULE, cds_results)
        assert cds_and_domain_names(hits) == [("cds1", "A"), ("cds2", "A"), ("cds2", "B")]


class TestBuildPredictions(unittest.TestCase):
    def setUp(self):
        self.rules = {THIS_RULE: build_rule(THIS_RULE), OTHER_RULE: build_rule(OTHER_RULE)}
        ruleset = SimpleNamespace(get_rule_by_name=self.rules.__getitem__)
        self.start_patch(patch.object(predictions, "get_signatures", return_value=SIGNATURES))
        self.get_ruleset = self.start_patch(patch.object(predictions, "get_ruleset",
                                                         return_value=ruleset))
        self.get_compound = self.start_patch(patch.object(predictions, "get_compound",
                                                          return_value=sentinel.compound))

    def start_patch(self, patcher):
        self.addCleanup(patcher.stop)
        return patcher.start()

    def build(self, cds_by_cluster, strictness="relaxed"):
        """ Builds predictions from rule detection results with the given
            CDS results for each protocluster
        """
        rule_results = SimpleNamespace(cds_by_cluster=cds_by_cluster)
        return build_predictions(rule_results, strictness)

    def test_prediction_details(self):
        protocluster = DummyProtocluster(start=10, end=500, product=THIS_RULE)
        cds_result = build_cds_result("cds1", hits=["A"], fired={THIS_RULE: {"A"}})

        predictions_built = self.build({protocluster: [cds_result]})

        assert len(predictions_built) == 1
        prediction = predictions_built[0]
        assert prediction.rule is self.rules[THIS_RULE]
        assert prediction.location == protocluster.location
        assert prediction.compound is sentinel.compound
        assert cds_and_domain_names(prediction.domain_hits) == [("cds1", "A")]
        self.get_compound.assert_called_once_with(THIS_RULE)

    def test_strictness_used(self):
        self.build({}, strictness="strict")
        self.get_ruleset.assert_called_once_with("strict")

    def test_one_per_protocluster(self):
        first = DummyProtocluster(start=10, end=500, product=THIS_RULE)
        second = DummyProtocluster(start=1000, end=1500, product=OTHER_RULE)

        predictions_built = self.build({
            first: [build_cds_result("cds1", hits=["A"], fired={THIS_RULE: {"A"}})],
            second: [build_cds_result("cds2", hits=["B"], fired={OTHER_RULE: {"B"}})],
        })

        # each prediction has the rule, location and hits of its own protocluster
        assert [prediction.rule_name for prediction in predictions_built] == [THIS_RULE,
                                                                              OTHER_RULE]
        assert [prediction.location for prediction in predictions_built] == [first.location,
                                                                             second.location]
        assert [cds_and_domain_names(prediction.domain_hits)
                for prediction in predictions_built] == [[("cds1", "A")], [("cds2", "B")]]

    def test_shared_cds(self):
        # a CDS in two overlapping protoclusters is shared by both, with an
        # entry for each rule, and each prediction only uses its own entry
        shared = build_cds_result("shared", hits=["A", "B"],
                                  fired={THIS_RULE: {"A"}, OTHER_RULE: {"B"}})
        first = DummyProtocluster(start=10, end=500, product=THIS_RULE)
        second = DummyProtocluster(start=300, end=800, product=OTHER_RULE)

        predictions_built = self.build({first: [shared], second: [shared]})

        assert [cds_and_domain_names(prediction.domain_hits)
                for prediction in predictions_built] == [[("shared", "A")], [("shared", "B")]]


class TestSubclusterPrediction(unittest.TestCase):
    def build(self, conditions="(A and B)", domain_hits=None):
        return SubclusterPrediction(rule=build_rule(THIS_RULE, conditions),
                                    location=sentinel.location, compound=sentinel.compound,
                                    domain_hits=domain_hits or [])

    def build_hit(self, cds_name, domain_name):
        return CDSDomainHit(domain_name=domain_name, domain_description=None,
                            domain_accession=None, cds_name=cds_name, evalue=1e-10, bitscore=50)

    def test_conditions_str(self):
        assert self.build("(A and B)").conditions_str == "A and B"
        # only the outermost brackets are removed
        assert self.build("((A or B) and C)").conditions_str == "(A or B) and C"
        assert self.build("A and B").conditions_str == "A and B"

    def test_grouping(self):
        cds1_a = self.build_hit("cds1", "A")
        cds2_a = self.build_hit("cds2", "A")
        cds1_b = self.build_hit("cds1", "B")
        prediction = self.build(domain_hits=[cds1_a, cds2_a, cds1_b])

        assert prediction.domain_hits_by_cds == {"cds1": [cds1_a, cds1_b], "cds2": [cds2_a]}
        assert prediction.domain_hits_by_domain == {"A": [cds1_a, cds2_a], "B": [cds1_b]}


class TestExtendsBeyond(unittest.TestCase):
    def extends(self, location, region=FeatureLocation(100, 500, 1)):
        prediction = SubclusterPrediction(rule=build_rule(THIS_RULE), location=location,
                                          compound=sentinel.compound, domain_hits=[])
        return prediction.extends_beyond(region)

    def test_contained(self):
        assert self.extends(FeatureLocation(150, 250, 1)) == (False, False)

    def test_same_as_region(self):
        # the end of a location is exclusive, so ending with the region doesn't extend it
        assert self.extends(FeatureLocation(100, 500, 1)) == (False, False)

    def test_extends_left(self):
        assert self.extends(FeatureLocation(50, 250, 1)) == (True, False)
        assert self.extends(FeatureLocation(99, 250, 1)) == (True, False)

    def test_extends_right(self):
        assert self.extends(FeatureLocation(150, 550, 1)) == (False, True)
        assert self.extends(FeatureLocation(150, 501, 1)) == (False, True)

    def test_extends_both(self):
        assert self.extends(FeatureLocation(50, 550, 1)) == (True, True)

    def test_cross_origin(self):
        # a region crossing the origin of a circular record of length 1000
        region = CompoundLocation([FeatureLocation(900, 1000, 1), FeatureLocation(0, 100, 1)])
        inside = CompoundLocation([FeatureLocation(950, 1000, 1), FeatureLocation(0, 50, 1)])
        assert self.extends(inside, region) == (False, False)
        left = CompoundLocation([FeatureLocation(850, 1000, 1), FeatureLocation(0, 50, 1)])
        assert self.extends(left, region) == (True, False)
        right = CompoundLocation([FeatureLocation(950, 1000, 1), FeatureLocation(0, 150, 1)])
        assert self.extends(right, region) == (False, True)
