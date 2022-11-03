# License: GNU Affero General Public License v3 or later
# A copy of GNU AGPL v3 should have been included in this software package in LICENSE.txt.

# for test files, silence irrelevant and noisy pylint warnings
# pylint: disable=use-implicit-booleaness-not-comparison,protected-access,missing-docstring

import unittest

from antismash.common.secmet import FeatureLocation
from antismash.modules.cluster_compare.data_structures import (
    Components,
    Hit,
    Mode,
    ReferenceArea,
    ReferenceCDS,
)


class DummyReferenceCDS(ReferenceCDS):
    counter = 0

    def __init__(self, name=None, function="other", components=None, location=None, start=None, strand=1):
        if name is None:
            DummyReferenceCDS.counter += 1
            name = f"test_ref_{DummyReferenceCDS.counter}"
        if components is None:
            components = {"secmet": [], "modules": []}
        if location is None:
            if start is None:
                start = 20
            location = FeatureLocation(start, start + 20, strand)
        super().__init__(name, function, components, location)


class DummyReferenceArea(ReferenceArea):
    def __init__(self, cds_mapping, cdses, accession="acc", start=5, end=500, products=None):
        if products is None:
            products = ["test"]
        super().__init__(accession, start, end, cds_mapping, cdses, products)


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
        components = Components({}, {}, {}, {})
        with self.assertRaises(AttributeError):
            components.something_invalid = "test"


class TestHit(unittest.TestCase):
    def test_restrictions(self):
        hit = Hit("acc", "id", None, 56., 250., 75., 1e-8)
        with self.assertRaises(AttributeError):
            hit.something_invalid = "test"
