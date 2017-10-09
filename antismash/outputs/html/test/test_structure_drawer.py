# License: GNU Affero General Public License v3 or later
# A copy of GNU AGPL v3 should have been included in this software package in LICENSE.txt.

# for test files, silence irrelevant and noisy pylint warnings
# pylint: disable=no-self-use,protected-access,missing-docstring

import os
import unittest

from helperlibs.wrappers.io import TemporaryDirectory

from antismash.outputs.html import structure_drawer


class TestLoadSmiles(unittest.TestCase):
    def test_load(self):
        smiles = structure_drawer.load_smiles()
        assert len(smiles) == 151
        # test a random selection
        assert smiles['23DHB'] == 'Oc1c(O)cccc1C(=O)'
        # and make sure none are empty
        for key, smile in smiles.items():
            assert smile, "smile %s has no smile" % key


class TestImageGeneration(unittest.TestCase):
    def test_depict_ectoine(self):
        with TemporaryDirectory(change=True) as temp:
            assert structure_drawer.depict_smile(0, "CC1=NCCC(N1)C(=O)O", temp)
            assert os.path.exists("genecluster0.smi")
            assert os.path.exists("genecluster0.png")
            assert os.path.exists("genecluster0_icon.png")
