# License: GNU Affero General Public License v3 or later
# A copy of GNU AGPL v3 should have been included in this software package in LICENSE.txt.

# for test files, silence irrelevant and noisy pylint warnings
# pylint: disable=no-self-use,protected-access,missing-docstring

import os
import unittest

from helperlibs.wrappers.io import TemporaryDirectory

from antismash.modules.nrps_pks import structure_drawer


class TestLoadSmiles(unittest.TestCase):
    def test_load(self):
        smiles = structure_drawer.load_smiles()
        assert len(smiles) == 152
        # test a random selection
        assert smiles['23DHB'] == 'Oc1c(O)cccc1C(=O)'
        # and make sure none are empty
        for key, smile in smiles.items():
            assert smile, "smile %s has no smile" % key
