# License: GNU Affero General Public License v3 or later
# A copy of GNU AGPL v3 should have been included in this software package in LICENSE.txt.

# for test files, silence irrelevant and noisy pylint warnings
# pylint: disable=use-implicit-booleaness-not-comparison,protected-access,missing-docstring

import unittest
from antismash.modules.nrps_pks.smiles_generator import (
    Atom,
    Bonds,
    gen_smiles_from_pksnrps as gen_smiles,
    get_all_smiles,
    methylate,
)


class TestGenerator(unittest.TestCase):
    components = {
        "mal": ("mal", "mal", []),
        "ohmal": ("mal", "ohmal", ["PKS_KR"]),
        "X": ("X", "X", []),
    }

    def test_empty(self):
        assert gen_smiles([]) == ""

    def test_single_nrp(self):
        assert gen_smiles([self.components["X"]])

    def test_single_ala(self):
        assert gen_smiles([("ala", "ala", [])]) == "NC(C)C(=O)O"
        # and a methylated variant
        assert gen_smiles([("ala", "me-ala", ["cMT"])]) == "NC(C(C))C(=O)O"

    def test_special_cases(self):
        # ending with a mal variant appends a pks-end2
        components = [self.components["ohmal"], self.components["mal"], ("mal", "ccmal", [])]
        smiles = gen_smiles(components)
        components.append(("pks-end2", "pks-end2", []))
        assert smiles == "CC(O)CC(=O)C=CC(=O)O" == gen_smiles(components)

    def test_mixed_mal(self):
        components = [self.components["mal"], self.components["X"], self.components["mal"]]
        assert gen_smiles(components) == "CC(=O)NC([*])C(=O)CC(=O)C(=O)O"


class TestAtom(unittest.TestCase):
    def test_simple_smiles(self):
        assert Atom("C").to_smiles() == "C"

    def test_bond_smiles(self):
        assert Atom("O", bonds_to_left=2).to_smiles() == "=O"
        assert Atom("C", bonds_to_left=3).to_smiles() == "#C"

    def test_identifier_smiles(self):
        assert Atom("C", identifier="2").to_smiles() == "C2"

    def test_branch_smiles(self):
        atom = Atom("C")
        atom.branches.append([Atom("C"), Atom("O")])
        atom.branches.append([Atom("N")])
        assert atom.to_smiles() == "C(CO)(N)"

    def test_combo(self):
        atom = Atom("C", identifier="1")
        atom.branches.append([Atom("C"), Atom("C", bonds_to_left=2)])
        atom.branches.append([Atom("O"), Atom("C")])
        assert atom.to_smiles() == "C1(C=C)(OC)"

    def test_nesting(self):
        atom = Atom("N")
        branching_atom = Atom("C")
        atom.branches.append([Atom("O"), branching_atom])
        branching_atom.branches.append([Atom("N"), Atom("C")])
        assert atom.to_smiles() == "N(OC(NC))"


class TestBonds(unittest.TestCase):
    def test_typical(self):
        smiles = "NC(C(=O)C)C(=O)O"
        bonds = Bonds(smiles)
        assert bonds.to_smiles() == smiles

    def test_multi_char_atoms(self):
        smiles = "CClC"
        bonds = Bonds(smiles)
        assert bonds.to_smiles() == smiles
        assert bonds._top_level_atoms[1].symbol == "Cl"

    def test_non_single_bonds(self):
        smiles = "C=CC#C"
        assert Bonds(smiles).to_smiles() == smiles

    def test_rings(self):
        smiles = "Cc1cc(O)cc(O)c1"
        bonds = Bonds(smiles)
        assert bonds.to_smiles() == smiles

    def test_multiple_references(self):
        smiles = "C1OC2OC12"
        bonds = Bonds(smiles)
        atoms = bonds._top_level_atoms
        assert atoms[0].references_in == [atoms[-1]]
        assert atoms[2].references_in == [atoms[-1]]
        assert atoms[-1].references_out == ["1", "2"]
        assert bonds.to_smiles() == smiles

    def test_iteration(self):
        bonds = Bonds("NC(C(=O)C)C(=O)O")
        assert [atom.symbol for atom in bonds] == ["N", "C", "C", "O", "C", "C", "O", "O"]

    def test_all_smiles(self):
        for key, val in get_all_smiles().items():
            assert Bonds(val).to_smiles() == val, key


def get_bond_counts(bonds):
    counts = []

    def per_atom_bonds(atom):
        bonds = [atom.available_bonds]
        if not atom.branches:
            return bonds
        children = []
        for branch in atom.branches:
            for child in branch:
                children.extend(per_atom_bonds(child))
        bonds.append(children)
        return bonds

    for atom in bonds._top_level_atoms:
        counts.extend(per_atom_bonds(atom))
    return counts


class TestBondCounts(unittest.TestCase):
    def test_dht(self):
        bonds = get_bond_counts(Bonds("NC(C(=O)C)C(=O)O"))
        #                N  C  (C (=O)  C)  C (=O)  O
        assert bonds == [2, 1, [0, [0], 3], 0, [0], 1]

    def test_double_bond_same_level(self):
        bonds = get_bond_counts(Bonds("C=C"))
        assert bonds == [2, 2]

    def test_ring(self):
        bonds = get_bond_counts(Bonds("C1NNC1"))
        #                C1 N  N  C 1
        assert bonds == [2, 1, 1, 2]

    def test_aromatic_carbon_rings(self):
        bonds = get_bond_counts(Bonds("c1cccc1"))
        assert bonds == [0, 0, 0, 0, 0]


class TestMethylation:
    def test_c(self):
        front = "NC(C(=O)C"
        back = ")C(=O)O"
        assert methylate(front + back, "C") == front + "(C)" + back

    def test_n(self):
        front = "NC(CCCN"
        back = ")C(=O)O"
        assert methylate(front + back, "N") == front + "(C)" + back

    def test_o(self):
        front = "NC(CCCC(=O)O"
        back = ")C(=O)O"
        assert methylate(front + back, "O") == front + "(C)" + back

    def test_ring_unmethylated(self):
        ring = "Nc1cccc1C(=O)O"
        assert methylate(ring, "C") == ring

    def test_last_o_skipped(self):
        smiles = "NC(CCCN)C(=O)O"
        assert methylate(smiles, "O") == smiles
