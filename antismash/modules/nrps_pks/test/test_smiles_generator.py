# License: GNU Affero General Public License v3 or later
# A copy of GNU AGPL v3 should have been included in this software package in LICENSE.txt.

# for test files, silence irrelevant and noisy pylint warnings
# pylint: disable=no-self-use,protected-access,missing-docstring

import unittest
from antismash.modules.nrps_pks.smiles_generator import gen_smiles_from_pksnrps as gen_smiles


class TestGenerator(unittest.TestCase):
    def test_empty(self):
        assert gen_smiles("") == ""

    def test_single_nrp(self):
        assert gen_smiles("(nrp)") == ""

    def test_single_ala(self):
        assert gen_smiles("(ala)") == ""

    def test_special_cases(self):
        # all mal variants and and end is ccmal appends a pks-end2
        polymer = "(ohmal - mal - ccmal)"
        assert gen_smiles(polymer) == "CC(O)CC(=O)C=CC(=O)(O)"
        assert gen_smiles("(ohmal - mal - ccmal - pks-end2)") == gen_smiles(polymer)

        # all mal variants and start is mal converts the first to pks-start1
        polymer = "(mal - ohmal - mal)"
        assert gen_smiles(polymer) == "CCC(O)CC(=O)"
        assert gen_smiles("(pks-start1 - ohmal - mal)") == gen_smiles(polymer)

        # pk in polymer and last is a mal variant removes the monomer after
        # the first pk and adds a pks-end1 at the end
        polymer = "(pk - nrp - pk - nrp - mal)"
        assert gen_smiles(polymer) == "C([*])C(-O)C([*])C(-O)NC([*])C(=O)CC(=O)C(C)C(=O)(O)"
        assert gen_smiles("(pk - pk - nrp - mal - pks-end1)") == gen_smiles(polymer)

    def test_mixed_mal(self):
        assert gen_smiles("(mal - nrp - mal)") == "CC(=O)NC([*])C(=O)CC(=O)"
