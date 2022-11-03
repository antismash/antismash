# License: GNU Affero General Public License v3 or later
# A copy of GNU AGPL v3 should have been included in this software package in LICENSE.txt.

# for test files, silence irrelevant and noisy pylint warnings
# pylint: disable=use-implicit-booleaness-not-comparison,protected-access,missing-docstring

import unittest

from antismash.common.secmet.features import PFAMDomain, FeatureLocation
from antismash.common.secmet.qualifiers import GOQualifier


class TestConversion(unittest.TestCase):
    def test_pfam_domain(self):
        original = PFAMDomain(FeatureLocation(2, 5), description="test",
                              protein_location=FeatureLocation(5, 10), identifier="PF00002.17",
                              domain="p450", tool="toolname", locus_tag="dummyCDS")
        original.domain_id = "domain_id"
        original.database = "db"
        original.detection = "someprogram"
        original.evalue = 1e-5
        original.score = 5.
        original.locus_tag = "locus"
        original.label = "somelabel"
        original.translation = "ARNDCQ"
        original.gene_ontologies = GOQualifier({'GO:0004871': 'signal transducer activity',
                                                'GO:0007165': 'signal transduction',
                                                'GO:0016020': 'membrane'})
        new = PFAMDomain.from_biopython(original.to_biopython()[0])
        for slot in ["tool", "domain_id", "database", "detection",
                     "evalue", "score", "locus_tag", "label", "translation", "domain",
                     "protein_location", "identifier", "version"]:
            assert getattr(original, slot) == getattr(new, slot)
        assert original.gene_ontologies.go_entries == new.gene_ontologies.go_entries
        assert original.full_identifier == new.full_identifier

    def test_bad_pfam_domain(self):
        protein_location = FeatureLocation(5, 10)
        with self.assertRaisesRegex(TypeError, "PFAMDomain description must be a string"):
            PFAMDomain(FeatureLocation(2, 5), description=None, protein_location=protein_location,
                       identifier="PF00002", tool="test", locus_tag="dummy")
        with self.assertRaisesRegex(TypeError, "Domain must be given domain as a string"):
            PFAMDomain(FeatureLocation(2, 5), description="desc", protein_location=protein_location,
                       identifier="PF00002", domain=5, tool="test", locus_tag="dummy")
        for ident in ["PF0002", "FAKE003", "PF", "PF000003", "PF00003.a"]:
            with self.assertRaisesRegex(ValueError, "invalid"):
                PFAMDomain(FeatureLocation(2, 5), description="desc", protein_location=protein_location,
                           identifier=ident, tool="test", locus_tag="dummy")
