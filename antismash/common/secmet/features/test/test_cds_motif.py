# License: GNU Affero General Public License v3 or later
# A copy of GNU AGPL v3 should have been included in this software package in LICENSE.txt.

# for test files, silence irrelevant and noisy pylint warnings
# pylint: disable=use-implicit-booleaness-not-comparison,protected-access,missing-docstring

import unittest

from Bio.SeqFeature import SeqFeature

from antismash.common.secmet.features import CDSMotif, FeatureLocation
from antismash.common.secmet.features.cds_motif import ExternalCDSMotif


class TestConversion(unittest.TestCase):
    def test_conversion(self):
        prot_loc = FeatureLocation(1, 2)
        original = CDSMotif(FeatureLocation(2, 5), tool="test", locus_tag="locus", protein_location=prot_loc)
        assert original.tool == "test"
        assert original.created_by_antismash
        assert original.locus_tag == "locus"
        assert original.protein_location == prot_loc

        bio_features = original.to_biopython()
        assert len(bio_features) == 1
        new = CDSMotif.from_biopython(bio_features[0])
        assert new.tool == original.tool == "test"
        assert new.locus_tag == original.locus_tag == "locus"
        assert new.protein_location == prot_loc
        assert new.created_by_antismash

    def test_non_antismash_motif(self):
        original = ExternalCDSMotif(FeatureLocation(7, 10), {})
        assert not original.created_by_antismash

        bio_features = original.to_biopython()
        assert len(bio_features) == 1, bio_features
        new = CDSMotif.from_biopython(bio_features[0])
        assert isinstance(new, ExternalCDSMotif)
        assert new.tool == original.tool
        assert not new.created_by_antismash

    def test_non_antismash_motif_from_raw(self):
        original = SeqFeature(FeatureLocation(7, 10))
        original.qualifiers["stuff"] = ["thing"]

        motif = CDSMotif.from_biopython(original)
        assert isinstance(motif, ExternalCDSMotif)
        assert motif.tool == "external"
        assert not motif.created_by_antismash
        assert motif.domain_id is None

        # add a domain_id so a Record can use the motif
        motif.domain_id = "testname"

        new = motif.to_biopython()[0]
        # generated domain_id should not be kept for non-antismash features
        assert "domain_id" not in new.qualifiers
        assert new.qualifiers == original.qualifiers
