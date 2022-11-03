# License: GNU Affero General Public License v3 or later
# A copy of GNU AGPL v3 should have been included in this software package in LICENSE.txt.

# for test files, silence irrelevant and noisy pylint warnings
# pylint: disable=use-implicit-booleaness-not-comparison,protected-access,missing-docstring

import unittest

from antismash.modules.lanthipeptides import rodeo


class TestRodeo(unittest.TestCase):
    def test_lanscout(self):
        assert rodeo.lanscout("ITSISLCTPGCKTGALMGCNMKTATCHCSIHVSK") == (7, [2, 1, 2, 2, 0])
