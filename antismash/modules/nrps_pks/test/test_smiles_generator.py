# License: GNU Affero General Public License v3 or later
# A copy of GNU AGPL v3 should have been included in this software package in LICENSE.txt.

# for test files, silence irrelevant and noisy pylint warnings
# pylint: disable=no-self-use,protected-access,missing-docstring

import unittest
from antismash.modules.nrps_pks.smiles_generator import (
    Atom,
    Bonds,
    gen_smiles_from_pksnrps as gen_smiles,
    load_smiles,
)


class TestGenerator(unittest.TestCase):
    def test_empty(self):
        assert gen_smiles("") == ""

    def test_single_nrp(self):
        assert gen_smiles("(X)")

    def test_single_ala(self):
        assert gen_smiles("(ala)")

    def test_special_cases(self):
        # all mal variants and and end is ccmal appends a pks-end2
        polymer = "(ohmal - mal - ccmal)"
        assert gen_smiles(polymer) == "CC(O)CC(=O)C=CC(=O)O"
        assert gen_smiles("(ohmal - mal - ccmal - pks-end2)") == gen_smiles(polymer)

        # all mal variants and start is mal converts the first to pks-start1
        polymer = "(mal - ohmal - mal)"
        assert gen_smiles(polymer) == "CCC(O)CC(=O)"
        assert gen_smiles("(pks-start1 - ohmal - mal)") == gen_smiles(polymer)

        # pk in polymer and last is a mal variant removes the monomer after
        # the first pk and adds a pks-end1 at the end
        polymer = "(pk - X - pk - X - mal)"
        assert gen_smiles(polymer) == "C([*])C(-O)C([*])C(-O)NC([*])C(=O)CC(=O)C(C)C(=O)O"
        assert gen_smiles("(pk - pk - X - mal - pks-end1)") == gen_smiles(polymer)

    def test_mixed_mal(self):
        assert gen_smiles("(mal - X - mal)") == "CC(=O)NC([*])C(=O)CC(=O)"


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
        for key, val in sorted(load_smiles().items()):
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
