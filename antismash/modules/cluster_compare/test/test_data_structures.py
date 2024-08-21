# License: GNU Affero General Public License v3 or later
# A copy of GNU AGPL v3 should have been included in this software package in LICENSE.txt.

# for test files, silence irrelevant and noisy pylint warnings
# pylint: disable=use-implicit-booleaness-not-comparison,protected-access,missing-docstring

import unittest
from unittest.mock import patch

from antismash.common import json
from antismash.common.secmet import FeatureLocation
from antismash.modules.cluster_compare.data_structures import (
    Components,
    Hit,
    Mode,
    ReferenceArea,
    ReferenceProtocluster,
    ReferenceRegion,
    ReferenceCDS,
)


class DummyReferenceCDS(ReferenceCDS):
    counter = 0

    def __init__(self, name=None, function="other", components=None, location=None, start=None, strand=1,
                 **kwargs):
        if name is None:
            DummyReferenceCDS.counter += 1
            name = f"test_ref_{DummyReferenceCDS.counter}"
        if components is None:
            components = {"secmet": [], "modules": []}
        if location is None:
            if start is None:
                start = 20
            location = FeatureLocation(start, start + 20, strand)
        super().__init__(name, function, components, location, **kwargs)


class DummyReferenceArea(ReferenceArea):
    def __init__(self, cds_mapping=None, cdses=None, accession="acc", start=5, end=500, products=None,
                 **kwargs,
                 ):
        if products is None:
            products = ["test"]
        if cdses is None:
            cdses = {name: DummyReferenceCDS() for name in cds_mapping.values()} \
                    if cds_mapping else {"A": DummyReferenceCDS()}
        if cds_mapping is None:
            cds_mapping = {str(i): cds.name for i, cds in enumerate(cdses.values())}

        super().__init__(accession, start, end, cds_mapping, cdses, products,
                         **kwargs)


class TestMode(unittest.TestCase):
    def test_restrictions(self):
        with self.assertRaises(AttributeError):
            Mode.something_that_doesnt_exist  # pylint: disable=no-member,pointless-statement

    def test_names(self):
        mode = Mode.REFERENCE_IN_QUERY
        assert mode.name == "REFERENCE_IN_QUERY"
        assert str(mode) == "RiQ"


class TestComponents(unittest.TestCase):
    def test_restrictions(self):
        components = Components({}, {}, {}, {}, {})
        with self.assertRaises(AttributeError):
            components.something_invalid = "test"


class TestJsonConversion(unittest.TestCase):
    def test_reference_area(self):
        area = DummyReferenceArea(start=350, end=450)
        expected = {
            "accession": area.accession,
            "start": area.start,
            "end": area.end,
            "draw_end": 450,
            "cdses": {key: cds.to_json() for key, cds in area.cdses.items()},
            "cds_mapping": area.cds_mapping,
            "products": area.products,
        }
        # to be sure that no objects exist in those mappings and such
        assert json.loads(json.dumps(expected)) == expected
        assert area.to_json() == expected
        # and a cross-origin area
        expected.update({
            "start": 450,
            "end": 50,
            "draw_end": 500,
        })
        area = DummyReferenceArea(start=450, end=50, draw_end=500, cdses=area.cdses, cds_mapping=area.cds_mapping)
        assert area.to_json() == expected

    def test_reference_region(self):
        proto = unittest.mock.Mock()
        proto.to_json = lambda: {"some key": "some value"}

        def check(region, expected_draw_end):
            assert region.draw_end == expected_draw_end
            with patch.object(ReferenceProtocluster, "from_json", return_value=proto):
                converted = ReferenceRegion.from_json(
                    region.accession, json.loads(json.dumps(region.to_json())),
                    cds_mapping=region.cds_mapping,
                )
            assert converted.to_json() == region.to_json()
            # and just to be sure that to_json isn't empty
            assert converted.draw_end == region.draw_end

        args = {
            "accession": "acc", "start": 10, "end": 500, "protoclusters": [proto],
            "cdses": {"A": DummyReferenceCDS()}, "products": ["p1", "p2"],
            "cds_mapping": {"1": "A"}, "description": "desc", "organism": "org",
        }
        check(ReferenceRegion(**args), expected_draw_end=500)
        args["draw_end"] = 550
        check(ReferenceRegion(**args), expected_draw_end=550)

    def test_cds_simple(self):
        cds = DummyReferenceCDS()
        converted = ReferenceCDS.from_json(cds.name, json.loads(json.dumps(cds.to_json())))
        assert cds.location.start == cds.draw_start == converted.draw_start
        assert cds.to_json() == converted.to_json()

    def test_cds_overrides(self):
        # mimic a cross- or post-origin CDS
        cds = DummyReferenceCDS(draw_start=600, draw_end=700)
        assert cds.draw_start != cds.location.start
        converted = ReferenceCDS.from_json(cds.name, json.loads(json.dumps(cds.to_json())))
        assert cds.draw_start == converted.draw_start
        assert cds.to_json() == converted.to_json()


class TestHit(unittest.TestCase):
    def test_restrictions(self):
        hit = Hit("acc", "id", None, 56., 250., 75., 1e-8)
        with self.assertRaises(AttributeError):
            hit.something_invalid = "test"
