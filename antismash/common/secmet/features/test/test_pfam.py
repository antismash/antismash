# License: GNU Affero General Public License v3 or later
# A copy of GNU AGPL v3 should have been included in this software package in LICENSE.txt.

# for test files, silence irrelevant and noisy pylint warnings
# pylint: disable=no-self-use,protected-access,missing-docstring

import unittest

from antismash.common.secmet.features import PFAMDomain, FeatureLocation
from antismash.common.secmet.qualifiers import GOQualifier


class TestConversion(unittest.TestCase):
    def test_pfam_domain(self):
        original = PFAMDomain(FeatureLocation(2, 5), description="test",
                              protein_start=5, protein_end=10, identifier="PF00002.17",
                              domain="p450", tool="toolname")
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
                     "protein_start", "protein_end", "identifier", "version"]:
            assert getattr(original, slot) == getattr(new, slot)
        assert original.gene_ontologies.go_entries == new.gene_ontologies.go_entries
        assert original.full_identifier == new.full_identifier

    def test_bad_pfam_domain(self):
        with self.assertRaisesRegex(TypeError, "PFAMDomain description must be a string"):
            PFAMDomain(FeatureLocation(2, 5), description=None, protein_start=5,
                       protein_end=10, identifier="PF00002", tool="test")
        with self.assertRaisesRegex(TypeError, "Domain must be given domain as a string"):
            PFAMDomain(FeatureLocation(2, 5), description="desc", protein_start=5,
                       protein_end=10, identifier="PF00002", domain=5, tool="test")
        with self.assertRaisesRegex(ValueError, "A PFAMDomain protein location cannot end before it starts"):
            PFAMDomain(FeatureLocation(2, 5), description="desc", protein_start=10,
                       protein_end=5, identifier="PF00002", tool="test")
        with self.assertRaisesRegex(ValueError, "invalid literal for int()"):
            PFAMDomain(FeatureLocation(2, 5), description="desc", protein_start=10,
                       protein_end="nope", identifier="PF00002", tool="test")
        for ident in ["PF0002", "FAKE003", "PF", "PF000003", "PF00003.a"]:
            with self.assertRaisesRegex(ValueError, "invalid"):
                PFAMDomain(FeatureLocation(2, 5), description="desc", protein_start=10,
                           protein_end="nope", identifier=ident, tool="test")
