# License: GNU Affero General Public License v3 or later
# A copy of GNU AGPL v3 should have been included in this software package in LICENSE.txt.

# for test files, silence irrelevant and noisy pylint warnings
# pylint: disable=use-implicit-booleaness-not-comparison,protected-access,missing-docstring

import unittest

from antismash.modules.cluster_compare import components
from antismash.modules.cluster_compare.data_structures import Components, Mode

from .test_data_structures import DummyReferenceArea, DummyReferenceCDS


class TestGatherRef(unittest.TestCase):
    def setUp(self):
        self.ref = DummyReferenceArea({"1": "A"}, {"A": DummyReferenceCDS()})
        assert self.ref._components is None

    def test_empty(self):
        result = components.gather_reference_components(self.ref)
        assert result is self.ref._components is self.ref.get_component_data()
        assert components.gather_reference_components(self.ref) is result
        assert result == components.Components({}, {}, {}, {})

    def test_secmet(self):
        self.ref.cdses["A"].components["secmet"] = ["name1", "name2"]
        result = components.gather_reference_components(self.ref)
        assert result is self.ref._components is self.ref.get_component_data()
        assert components.gather_reference_components(self.ref) is result
        assert result.secmet == {"name1": 1, "name2": 1}
        assert not any([result.functions, result.nrps, result.pks])

    def test_nrps(self):
        module1 = {
            "complete": True,
            "type": "nrps",
            "domains": ["Condensation", "AMP-binding", "PP-binding"],
        }
        module2 = {
            "complete": False,
            "type": "nrps",
            "domains": ["AMP-binding", "PP-binding"],
        }
        self.ref.cdses["A"].components["modules"] = [module1, module2]
        result = components.gather_reference_components(self.ref)
        assert result.nrps == {tuple(module1["domains"]): 1}
        assert not any([result.functions, result.secmet, result.pks])

    def test_pks(self):
        module1 = {
            "complete": True,
            "type": "pks",
            "domains": ["PKS_KS", "AT", "PP-binding"],
        }
        module2 = {
            "complete": False,
            "type": "pks",
            "domains": ["AT", "PP-binding"],
        }
        self.ref.cdses["A"].components["modules"] = [module1, module2]
        result = components.gather_reference_components(self.ref)
        assert result.pks == {tuple(module1["domains"]): 1}
        assert not any([result.functions, result.secmet, result.nrps])

    def test_function(self):
        self.ref.cdses["A"].function = "biosynthetic"
        result = components.gather_reference_components(self.ref)
        assert result.functions == {"biosynthetic": 1}
        assert not any([result.pks, result.secmet, result.nrps])


class TestChunkComparison(unittest.TestCase):
    def test_directionality_members(self):
        ref = {"A": 1, "X": 1}
        query = {"A": 1}
        result = components.compare_combos(ref, query, Mode.REFERENCE_IN_QUERY)
        assert result == 0.5
        result = components.compare_combos(ref, query, Mode.QUERY_IN_REFERENCE)
        assert result == 1.0

    def test_directionality_counts(self):
        ref = {"A": 2}
        query = {"A": 1}
        result = components.compare_combos(ref, query, Mode.REFERENCE_IN_QUERY)
        assert result == 0.5
        result = components.compare_combos(ref, query, Mode.QUERY_IN_REFERENCE)
        assert result == 1.0

    def test_directionality_empty(self):
        empty = {}
        full = {"A": 1}
        result = components.compare_combos(empty, full, Mode.REFERENCE_IN_QUERY)
        assert result == 1.
        result = components.compare_combos(empty, full, Mode.QUERY_IN_REFERENCE)
        assert result == 0.
        result = components.compare_combos(full, empty, Mode.REFERENCE_IN_QUERY)
        assert result == 0.
        result = components.compare_combos(full, empty, Mode.QUERY_IN_REFERENCE)
        assert result == 1.


class TestComplete(unittest.TestCase):
    def test_function_only_ref(self):
        cds = DummyReferenceCDS()
        cds.components["function"] = "test"
        ref = DummyReferenceArea({"1": "A"}, {"A": cds})
        query = Components({}, {}, {}, {"test": 1})
        result = components.calculate_component_score(query, ref, Mode.REFERENCE_IN_QUERY)
        assert result is None
        result = components.calculate_component_score(query, ref, Mode.QUERY_IN_REFERENCE)
        assert result == 0.

    def test_mixed(self):
        ref = Components(nrps={("A", "B"): 1, ("B", "C"): 1},
                         pks={("E", "F"): 1, ("F",): 1},
                         secmet={"X": 1},
                         functions={"biosynthetic": 3, "regulatory": 1})
        query = Components(nrps={("A", "B"): 1},
                           pks={("E", "F"): 1, ("F",): 1},
                           secmet={"X": 1, "Y": 1},
                           functions={"biosynthetic": 2, "regulatory": 1, "transport": 2})
        result = components.compare(ref, query, Mode.REFERENCE_IN_QUERY)
        self.assertAlmostEqual(result, 0.833, places=3)
        result = components.compare(ref, query, Mode.QUERY_IN_REFERENCE)
        self.assertAlmostEqual(result, 0.700, places=3)
