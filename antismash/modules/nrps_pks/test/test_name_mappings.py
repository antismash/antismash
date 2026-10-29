# License: GNU Affero General Public License v3 or later
# A copy of GNU AGPL v3 should have been included in this software package in LICENSE.txt.

# for test files, silence irrelevant and noisy pylint warnings
# pylint: disable=use-implicit-booleaness-not-comparison,protected-access,missing-docstring

from collections import defaultdict
from operator import sub
import unittest

from antismash.modules.nrps_pks import name_mappings


class TestNameMappings(unittest.TestCase):
    def setUp(self) -> None:
        self.known = name_mappings.KNOWN_SUBSTRATES

    def test_long_name_unique(self) -> None:
        long_names = set()
        counts: dict[str, int] = defaultdict(int)

        for substrate in self.known:
            long_names.add(substrate.long)
            counts[substrate.short] += 1

        assert len(self.known) == len(long_names), counts

    def test_short_name_unique(self) -> None:
        short_names = set()
        counts: dict[str, int] = defaultdict(int)

        for substrate in self.known:
            short_names.add(substrate.short)
            counts[substrate.short] += 1

        assert len(self.known) == len(short_names), [(name, count) for name, count in counts.items() if count > 1]

    def test_loaded_substrate_metadata(self) -> None:
        substrate = name_mappings.get_substrate_by_name("alanine")
        assert (substrate.long, substrate.short, substrate.norine) == ("alanine", "Ala", "Ala")
        assert substrate.smiles == "NC(C)C(=O)O"
        assert substrate.stereochemical_smiles == "N[C@@H](C)C(=O)O"
        assert substrate.substrate_type == "alpha_amino_acid"
        assert substrate.parent is name_mappings.CATEGORIES_BY_ID["tiny_alpha_amino_acids"]
        assert name_mappings.SUBSTRATES_BY_ID["ala"] is substrate
        assert substrate.parent is not None
        assert any(child is substrate for child in substrate.parent.children)
        for name in ("ALANINE", "ala", "Ala"):
            assert name_mappings.get_substrate_by_name(name) is substrate

    def test_taxonomy_complete(self) -> None:
        ids = set()
        for substrate in self.known:
            assert substrate.id
            assert substrate.id not in ids
            ids.add(substrate.id)
            assert name_mappings.get_substrate_by_name(substrate.id) is substrate
            assert substrate.smiles
            assert substrate.substrate_type
            assert [ancestor.level for ancestor in substrate.ancestors] == [
                "small_cluster", "large_cluster", "gross_class",
            ]
            assert substrate.parent is not None
            assert any(child is substrate for child in substrate.parent.children)

    def test_pks_fragments(self) -> None:
        substrate = name_mappings.get_substrate_by_name("mal")
        assert substrate.smiles == "CC(=O)O"
        assert substrate.substrate_type == "polyketide_fragment"
        assert substrate.parent is name_mappings.CATEGORIES_BY_ID["carbonyl_polyketide_units"]
        # The prediction labels do not determine hydroxyl configuration or alkene geometry.
        for name in ("ohmal", "ccmal", "me-redmal"):
            assert name_mappings.get_substrate_by_name(name).stereochemical_smiles is None
        substrate = name_mappings.get_substrate_by_name("redmal")
        assert substrate.smiles == substrate.stereochemical_smiles == "CC"

    def test_unknown_substrate_structures(self) -> None:
        substrate = name_mappings.get_substrate_by_name("nrp")
        assert substrate.smiles == "NC([*])C(=O)O"
        assert substrate.stereochemical_smiles is None
        assert substrate.ancestors[-1].id == "unclassified"

        substrate = name_mappings.get_substrate_by_name("pk")
        assert substrate.smiles == "C([*])C(=O)O"
        assert substrate.stereochemical_smiles is None
        assert substrate.ancestors[-1].id == "unclassified"

    def test_substrate_serialisation(self) -> None:
        for substrate in self.known:
            data = substrate.to_json()
            restored = name_mappings.Substrate.from_json(data)
            assert restored.to_json() == data
            assert restored.parent is substrate.parent

    def test_category_serialisation(self) -> None:
        restored_categories = {}
        for category in name_mappings.KNOWN_CATEGORIES:
            data = category.to_json()
            assert "parent" not in data
            assert data["parent_id"] == (category.parent.id if category.parent else None)
            restored = name_mappings.Category.from_json(restored_categories, **data)
            assert restored == category
            assert restored.to_json() == data
            if category.parent is None:
                assert restored.parent is None
            else:
                assert restored.parent is restored_categories[category.parent.id]
            restored_categories[restored.id] = restored

    def test_smiles(self) -> None:
        for substrate in self.known:
            if substrate.smiles:
                if substrate.substrate_type.endswith("_amino_acid"):
                    assert substrate.smiles.startswith("N")
                if substrate.substrate_type.endswith("acid"):
                    # ignore a couple of weird cases
                    if substrate.short in ("nHcp", "Ahp"):
                        continue
                    # Preserve carboxylates in zwitterionic substrates.
                    assert substrate.smiles.endswith(("C(=O)O", "C(=O)[O-]")), substrate.short
