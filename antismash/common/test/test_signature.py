# License: GNU Affero General Public License v3 or later
# A copy of GNU AGPL v3 should have been included in this software package in LICENSE.txt.

# for test files, silence irrelevant and noisy pylint warnings
# pylint: disable=use-implicit-booleaness-not-comparison,protected-access,missing-docstring,consider-using-with

from unittest import TestCase

from antismash.common.signature import Signature


class TestSignature(TestCase):
    def test_whitespace_in_identifiers(self):
        for name in [" leading", "trailing ", " both "]:
            with self.assertRaisesRegex(ValueError, "cannot have leading or trailing"):
                Signature(name, "type", "description", 50, "dummy_path")
        assert Signature("good", "type", "description", 50, "dummy_path").name == "good"
