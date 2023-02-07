# License: GNU Affero General Public License v3 or later
# A copy of GNU AGPL v3 should have been included in this software package in LICENSE.txt.

# for test files, silence irrelevant and noisy pylint warnings
# pylint: disable=use-implicit-booleaness-not-comparison,protected-access,missing-docstring,too-many-public-methods

import json
import unittest

from antismash.common.secmet.features.protocluster import SideloadedProtocluster
from antismash.common.secmet.features.subregion import SideloadedSubRegion
from antismash.common.test.helpers import DummyRecord
from antismash.detection.sideloader import data_structures as structures


def dummy_tool():
    return structures.Tool("name", "vers", "desc", configuration={})


class TestSub(unittest.TestCase):
    def test_bad_location(self):
        with self.assertRaisesRegex(ValueError, "end must be greater"):
            structures.SubRegionAnnotation(10, 5, "label", dummy_tool(), {})

    def test_attributes(self):
        sub = structures.SubRegionAnnotation(5, 50, "label", dummy_tool(),
                                             {"a": ["first", "second"], "b": ["third"]})
        assert len(sub) == sub.end - sub.start == 45

    def test_json_conversion(self):
        sub = structures.SubRegionAnnotation(5, 50, "label", dummy_tool(), {"a": ["first"], "b": ["and", "second"]})
        as_json = json.loads(json.dumps(sub.to_json()))  # string conversion to match real use
        regenerated = structures.SubRegionAnnotation.from_json(as_json)
        assert sub == regenerated, f"{vars(sub)} != {vars(regenerated)}"

        as_json["details"] = {"a": "first", "b": ["and", "second"]}
        regenerated = structures.SubRegionAnnotation.from_json(as_json)
        assert regenerated.details["a"] == ["first"]  # single strings should be listified

        as_json.pop("tool")
        regenerated = structures.SubRegionAnnotation.from_schema_json(as_json, dummy_tool())
        assert sub == regenerated, f"{vars(sub)} != {vars(regenerated)}"

    def test_secmet_conversion(self):
        source = structures.SubRegionAnnotation(5, 50, "label", dummy_tool(), {"a": ["first"], "b": ["second"]})
        dest = source.to_secmet()
        assert isinstance(dest, SideloadedSubRegion)
        assert dest.label == "label"
        assert dest.tool == source.tool.name
        assert dest.location.start == 5
        assert dest.location.end == 50

    def test_equality(self):
        first = structures.SubRegionAnnotation(5, 50, "label", dummy_tool(), {"a": ["first"], "b": ["second"]})
        second = structures.SubRegionAnnotation(5, 50, "label", dummy_tool(), {"a": ["first"], "b": ["second"]})
        assert first == second
        assert first != second.details
        second.details["b"][0] = "econd"
        assert first != second


class TestCluster(unittest.TestCase):
    def test_bad_location(self):
        with self.assertRaisesRegex(ValueError, "end must be greater"):
            structures.ProtoclusterAnnotation(10, 5, "product", dummy_tool(), {})

        with self.assertRaisesRegex(ValueError, "absolute distance"):
            structures.ProtoclusterAnnotation(5, 10, "product", dummy_tool(), {}, neighbourhood_left=-1)
        with self.assertRaisesRegex(ValueError, "absolute distance"):
            structures.ProtoclusterAnnotation(5, 10, "product", dummy_tool(), {}, neighbourhood_right=-1)

        with self.assertRaisesRegex(ValueError, "cannot extend out of a record"):
            structures.ProtoclusterAnnotation(5, 10, "product", dummy_tool(), {}, neighbourhood_left=7)

    def test_attributes(self):
        cluster = structures.ProtoclusterAnnotation(25, 50, "product", dummy_tool(),
                                                    {"a": ["first"], "b": ["second"]},
                                                    neighbourhood_left=10, neighbourhood_right=5)
        assert cluster.start == 15
        cluster.neighbourhood_left += 3
        assert cluster.start == 12

        assert cluster.end == 55
        cluster.neighbourhood_right += 3
        assert cluster.end == 58

        assert len(cluster) == cluster.end - cluster.start

        assert cluster.label == cluster.product == "product"
        cluster.product = "something else"
        assert cluster.label == cluster.product == "something else"

        assert isinstance(cluster.tool, structures.Tool)
        assert cluster.tool.name == dummy_tool().name

    def test_json_conversion(self):
        cluster = structures.ProtoclusterAnnotation(5, 50, "product", dummy_tool(),
                                                    {"a": ["first"], "b": ["and", "second"]})
        as_json = json.loads(json.dumps(cluster.to_json()))  # string conversion to match real use
        regenerated = structures.ProtoclusterAnnotation.from_json(as_json)
        assert cluster == regenerated

        as_json["details"] = {"a": "first", "b": ["and", "second"]}
        regenerated = structures.ProtoclusterAnnotation.from_json(as_json)
        assert regenerated.details["a"] == ["first"]  # single strings should be listified

        as_json.pop("tool")
        regenerated = structures.ProtoclusterAnnotation.from_schema_json(as_json, dummy_tool())
        assert cluster == regenerated

    def test_secmet_conversion(self):
        source = structures.ProtoclusterAnnotation(25, 50, "product", dummy_tool(),
                                                   {"a": ["first"], "b": ["second"]},
                                                   neighbourhood_left=10, neighbourhood_right=5)
        dest = source.to_secmet()
        assert isinstance(dest, SideloadedProtocluster)
        assert dest.product == "product"
        assert dest.tool == source.tool.name
        assert dest.core_location.start == 25
        assert dest.core_location.end == 50
        assert dest.location.start == 15
        assert dest.location.end == 55
        assert dest.neighbourhood_range == 10

    def test_equality(self):
        first = structures.ProtoclusterAnnotation(5, 50, "product", dummy_tool(), {"a": ["first"]})
        second = structures.ProtoclusterAnnotation(5, 50, "product", dummy_tool(), {"a": ["first"]})
        assert first == second
        assert first != second.details
        second.details["a"][0] = "irst"
        assert first != second


class TestTool(unittest.TestCase):
    def test_simple(self):
        tool = structures.Tool("name", "2.5.6", "desc of tool", {"configA": ["valA1", "valA2"]})
        assert tool.name == "name"
        assert tool.version == "2.5.6"
        assert tool.description == "desc of tool"
        assert tool.configuration == {"configA": ["valA1", "valA2"]}

    def test_equality(self):
        first = structures.Tool("name", "2.5.6", "desc of tool", {"configA": ["valA1", "valA2"]})
        second = structures.Tool("name", "2.5.6", "desc of tool", {"configA": ["valA1", "valA2"]})
        assert first == second
        assert first != second.configuration
        second.configuration["configA"][1] = "valB2"
        assert first != second

    def test_json_conversion(self):
        source = structures.Tool("name", "2.5.6", "desc of tool", {"configA": ["valA1", "valA2"]})
        as_json = json.loads(json.dumps(source.to_json()))
        dest = structures.Tool.from_json(as_json)
        assert source.name == dest.name
        assert source.version == dest.version
        assert source.description == dest.description
        assert source.configuration == dest.configuration

    def test_names(self):
        for name in ["a_b", "a b", "a-b"]:
            assert structures.Tool(name, "version", "desc", {})
        for name in ["a/b", "a!", "a#", "$a", "a\tb"]:
            with self.assertRaisesRegex(ValueError, "invalid character"):
                structures.Tool(name, "version", "desc", {})


class TestResults(unittest.TestCase):
    def setUp(self):
        self.sub = structures.SubRegionAnnotation(5, 50, "label", dummy_tool(),
                                                  {"a": ["first"], "b": ["and", "second"]})
        self.cluster = structures.ProtoclusterAnnotation(5, 40, "product", dummy_tool(), {"a": ["first"]})
        self.results = structures.SideloadedResults("rec_name", [self.sub], [self.cluster])

    def test_regeneration(self):
        record = DummyRecord()
        record.id = "rec_name"
        as_json = json.loads(json.dumps(self.results.to_json()))
        regenerated = structures.SideloadedResults.from_json(as_json, record)
        assert regenerated.protoclusters == self.results.protoclusters
        assert regenerated.subregions == self.results.subregions
        assert regenerated.record_id == self.results.record_id

        as_json["schema_version"] = -1
        with self.assertRaisesRegex(ValueError, "Detection results have changed"):
            structures.SideloadedResults.from_json(as_json, record)

    def test_areas(self):
        assert self.results.get_areas() == [self.sub, self.cluster]
        self.sub.start = 6
        assert self.results.get_areas() == [self.cluster, self.sub]
        extra = structures.ProtoclusterAnnotation(2, 50, "product", dummy_tool(), {})
        self.results.protoclusters.append(extra)
        assert self.results.get_areas() == [extra, self.cluster, self.sub]

    def test_secmet(self):
        subs = self.results.get_predicted_subregions()
        assert len(subs) == 1
        assert isinstance(subs[0], SideloadedSubRegion)

        clusters = self.results.get_predicted_protoclusters()
        assert len(clusters) == 1
        assert isinstance(clusters[0], SideloadedProtocluster)
