# License: GNU Affero General Public License v3 or later
# A copy of GNU AGPL v3 should have been included in this software package in LICENSE.txt.

# for test files, silence irrelevant and noisy pylint warnings
# pylint: disable=no-self-use,protected-access,missing-docstring

import unittest

from antismash.common.secmet.features import CDSMotif, FeatureLocation


class TestConversion(unittest.TestCase):
    def test_tool_conversion(self):
        original = CDSMotif(FeatureLocation(2, 5), tool="test")
        assert original.tool == "test"
        assert original.created_by_antismash

        bio_features = original.to_biopython()
        assert len(bio_features) == 1
        new = CDSMotif.from_biopython(bio_features[0])
        assert new.tool == original.tool == "test"
        assert new.created_by_antismash

    def test_non_antismash_motif(self):
        original = CDSMotif(FeatureLocation(7, 10))
        assert original.tool is None
        assert not original.created_by_antismash

        bio_features = original.to_biopython()
        assert len(bio_features) == 1
        new = CDSMotif.from_biopython(bio_features[0])
        assert new.tool is None
        assert not new.created_by_antismash
