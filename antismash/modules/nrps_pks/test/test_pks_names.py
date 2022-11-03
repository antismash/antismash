# License: GNU Affero General Public License v3 or later
# A copy of GNU AGPL v3 should have been included in this software package in LICENSE.txt.

# for test files, silence irrelevant and noisy pylint warnings
# pylint: disable=use-implicit-booleaness-not-comparison,protected-access,missing-docstring

import unittest
from antismash.modules.nrps_pks.pks_names import get_long_form, get_short_form


class TestConversion(unittest.TestCase):
    def test_known_short(self):
        assert get_long_form("mal") == "Malonyl-CoA"
        assert get_long_form("isobut") == "Isobutyryl-CoA"

    def test_unknown_short(self):
        assert get_long_form("foo") == "foo"
        assert get_long_form("foo", "other") == "other"

    def test_known_long(self):
        assert get_short_form("Malonyl-CoA") == "mal"
        assert get_short_form("Isobutyryl-CoA") == "isobut"

    def test_unknown_long(self):
        assert get_short_form("unknown") == "unknown"
        assert get_short_form("unknown", "other") == "other"
