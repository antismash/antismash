# License: GNU Affero General Public License v3 or later
# A copy of GNU AGPL v3 should have been included in this software package in LICENSE.txt.

# for test files, silence irrelevant and noisy pylint warnings
# pylint: disable=use-implicit-booleaness-not-comparison,protected-access,missing-docstring

import unittest
from unittest.mock import patch

from antismash.common.layers import RecordLayer, RegionLayer
from antismash.common.secmet.locations import FeatureLocation, CompoundLocation
from antismash.common.test.helpers import DummyRecord, DummyRegion
from antismash.config import build_config, destroy_config


@patch.object(RegionLayer, "find_plugins_for_region")
class TestRegionDescription(unittest.TestCase):
    def setUp(self):
        record = DummyRecord(seq="A" * 200)
        self.region = DummyRegion()
        record.add_region(self.region)
        with patch.object(RegionLayer, "find_plugins_for_region"):
            self.record_layer = RecordLayer(record, None, build_config([]))

    def tearDown(self):
        destroy_config()

    def test_coordinates(self, _patched_find):
        self.region.location = FeatureLocation(5, 150, 1)
        region_layer = RegionLayer(self.record_layer, self.region)
        assert "Location: 6 - 150 nt" in region_layer.description_text()  # 1-indexed

    def test_coordinates_cross_origin(self, _patched_find):
        self.region.location = CompoundLocation([
            FeatureLocation(150, 200, 1),
            FeatureLocation(0, 5, 1),
        ])
        region_layer = RegionLayer(self.record_layer, self.region)
        assert "Location: 151 - 5 nt" in region_layer.description_text()  # 1-indexed
