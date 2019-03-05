# License: GNU Affero General Public License v3 or later
# A copy of GNU AGPL v3 should have been included in this software package in LICENSE.txt.

# for test files, silence irrelevant and noisy pylint warnings
# pylint: disable=no-self-use,protected-access,missing-docstring

import unittest

from Bio.SeqFeature import SeqFeature

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

    def test_non_antismash_motif_from_raw(self):
        original = SeqFeature(FeatureLocation(7, 10))
        original.qualifiers["stuff"] = ["thing"]

        motif = CDSMotif.from_biopython(original)
        assert motif.tool is None
        assert not motif.created_by_antismash
        assert motif.domain_id is None

        # add a domain_id so a Record can use the motif
        motif.domain_id = "testname"

        new = motif.to_biopython()[0]
        # generated domain_id should not be kept for non-antismash features
        assert "domain_id" not in new.qualifiers
        assert new.qualifiers == original.qualifiers
