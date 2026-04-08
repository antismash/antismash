# License: GNU Affero General Public License v3 or later
# A copy of GNU AGPL v3 should have been included in this software package in LICENSE.txt.

# for test files, silence irrelevant and noisy pylint warnings
# pylint: disable=use-implicit-booleaness-not-comparison,protected-access,missing-docstring,too-many-public-methods

import unittest
from unittest.mock import patch

from antismash.common.secmet.test.helpers import (
    DummyRecord,
    DummyRegion,
)
from antismash.common.html_renderer import FileTemplate as _FileTemplate
from antismash.common.layers import RegionLayer
from antismash.detection.sideloader.data_structures import (
    SideloadedResults,
    SubRegionAnnotation,
    Tool,
)
from antismash.detection.sideloader import html_output


def create_dummy_tool():
    return Tool("name", "version", "description", {})


class TestAreaFilter(unittest.TestCase):
    def check(self, region, results, expected_areas):
        class MockTemplate(_FileTemplate):
            def render(self, areas=None, **kwargs):
                assert {area.label for area in areas} == {area.label for area in expected_areas}

        with patch.object(html_output, "FileTemplate", MockTemplate):
            html_output.generate_html(region, results, None, None)

    def test_linear(self):
        all_areas = [
            SubRegionAnnotation(10, 20, "label1", create_dummy_tool(), {}),
            SubRegionAnnotation(35, 45, "label2", create_dummy_tool(), {}),
            SubRegionAnnotation(50, 55, "label3", create_dummy_tool(), {}),
        ]
        results = SideloadedResults("record id", subregions=all_areas, protoclusters=[])

        region = DummyRegion(start=30, end=60, subregions=[])
        record = DummyRecord(length=100, features=[region], circular=False)
        with patch.object(RegionLayer, "find_plugins_for_region", return_value=[]):
            self.check(RegionLayer(record, region), results, expected_areas=all_areas[1:])

    def test_circular(self):
        length = 100
        all_areas = [
            SubRegionAnnotation(80, 10, "cross-origin", create_dummy_tool(), {},
                                circular_origin=length),
            SubRegionAnnotation(15, 20, "early", create_dummy_tool(), {}),
            SubRegionAnnotation(25, 30, "mid", create_dummy_tool(), {}),
            SubRegionAnnotation(50, 60, "late", create_dummy_tool(), {}),
        ]
        results = SideloadedResults("record id", subregions=all_areas, protoclusters=[])

        region_mid = DummyRegion(start=25, end=40, subregions=[])
        region_origin = DummyRegion(start=45, end=20, circular_wrap_point=length, subregions=[])
        record = DummyRecord(length=length, features=[region_mid, region_origin], circular=True)

        with patch.object(RegionLayer, "find_plugins_for_region", return_value=[]):
            self.check(RegionLayer(record, region_mid), results, expected_areas=[all_areas[2]])
        with patch.object(RegionLayer, "find_plugins_for_region", return_value=[]):
            self.check(RegionLayer(record, region_origin), results, expected_areas=[
                all_areas[0], all_areas[1], all_areas[-1],
            ])
