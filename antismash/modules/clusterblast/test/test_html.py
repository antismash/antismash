# License: GNU Affero General Public License v3 or later
# A copy of GNU AGPL v3 should have been included in this software package in LICENSE.txt.

# for test files, silence irrelevant and noisy pylint warnings
# pylint: disable=use-implicit-booleaness-not-comparison,protected-access,missing-docstring

import unittest
from unittest.mock import Mock, patch

from antismash.common.layers import OptionsLayer, RecordLayer, RegionLayer
from antismash.common.secmet.locations import CompoundLocation, FeatureLocation
from antismash.common.secmet.test.helpers import DummyCDS, DummyRecord, DummyRegion
from antismash.common.test.helpers import get_simple_options
from antismash.modules import clusterblast


class TestHTMLReuse(unittest.TestCase):
    """ Catches the case where only the options were used to determine
        whether HTML sections should be generated
    """
    def setUp(self):
        self.options = OptionsLayer(get_simple_options(clusterblast, []), modules=[clusterblast])
        region = DummyRegion()
        self.record = RecordLayer(DummyRecord(features=[region]), results=None, options=self.options)
        self.region = RegionLayer(self.record, region)

        self.results = Mock(spec=["general", "knowncluster", "subcluster"],
                            general=None, knowncluster=None, subcluster=None)

    def build_result(self):
        result = clusterblast.results.GeneralResults(self.record.id, data_version=1)
        result.region_results = {self.region.get_region_number() - 1: Mock(ranking=[])}
        return result

    def check(self, *, expected_call_count=1):
        with patch.object(clusterblast.html_output, "generate_div", return_value="") as patched:
            sections = clusterblast.generate_html(self.region, self.results, self.record, self.options)
            assert len(patched.mock_calls) == expected_call_count
        return sections

    def test_general(self):
        self.results.general = self.build_result()
        sections = self.check()
        assert len(sections.detail_sections) == 1
        assert sections.detail_sections[0].label == "ClusterBlast"

    def test_known(self):
        self.results.knowncluster = self.build_result()
        sections = self.check()
        assert len(sections.detail_sections) == 1
        assert sections.detail_sections[0].label == "KnownClusterBlast"

    def test_sub(self):
        self.results.subcluster = self.build_result()
        sections = self.check()
        assert len(sections.detail_sections) == 1
        assert sections.detail_sections[0].label == "SubClusterBlast"

    def test_all(self):
        self.results.knowncluster = self.build_result()
        self.results.general = self.build_result()
        self.results.subcluster = self.build_result()
        sections = self.check(expected_call_count=3)
        assert len(sections.detail_sections) == 3
        assert {section.label for section in sections.detail_sections} == {
            "ClusterBlast",
            "KnownClusterBlast",
            "SubClusterBlast",
        }


class TestQueryJSON(unittest.TestCase):
    def test_full_record_with_cross_origin(self):
        region = DummyRegion(start=0, end=100)
        standard = DummyCDS(start=50, end=80)
        cross_origin = DummyCDS(location=CompoundLocation([
            FeatureLocation(90, 100, 1),
            FeatureLocation(0, 20, 1),
        ]))
        region.add_cds(standard)
        region.add_cds(cross_origin)
        result = clusterblast.html_output.QueryJSON.from_region(region)
        assert len(result.genes) == 2
        for initial, built in zip([cross_origin, standard], result.genes):
            assert initial.get_name() == built.locus_tag
            assert initial.start == built.start
        assert result.genes[0].end == 120  # since it's pushed over the origin
