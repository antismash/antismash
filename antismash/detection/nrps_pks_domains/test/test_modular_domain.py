# License: GNU Affero General Public License v3 or later
# A copy of GNU AGPL v3 should have been included in this software package in LICENSE.txt.

# for test files, silence irrelevant and noisy pylint warnings
# pylint: disable=use-implicit-booleaness-not-comparison,protected-access,missing-docstring

import unittest

from antismash.common.secmet.features import (
    AntismashDomain,
    FeatureLocation,
)

from antismash.detection.nrps_pks_domains.modular_domain import ModularDomain, TOOL


class TestConversion(unittest.TestCase):
    def setUp(self):
        self.protein_location = FeatureLocation(0, 1)
        self.domain = ModularDomain(FeatureLocation(1, 3, 1), locus_tag="locus",
                                    protein_location=self.protein_location)
        self.domain.domain_id = "test1"

    def test_conversion(self):
        domain = self.domain
        domain.subtypes = ["subtest", "more"]
        domain.specificity = ["a", "c", "f"]
        domain.asf.add("first")
        domain.asf.add("second")
        assert domain.tool == TOOL
        assert domain.created_by_antismash
        assert domain.locus_tag == "locus"

        bio = domain.to_biopython()
        assert len(bio) == 1
        assert bio[0].qualifiers["aSTool"] == ["nrps_pks_domains"]
        assert bio[0].qualifiers["tool"] == ["antismash"]
        new_domain = ModularDomain.from_biopython(bio[0])
        assert new_domain.subtypes == domain.subtypes == ["subtest", "more"]
        assert new_domain.specificity == domain.specificity == ["a", "c", "f"]
        assert new_domain.asf.hits == domain.asf.hits
        assert new_domain.asf.hits == ["first", "second"]
        assert new_domain.tool == domain.tool == TOOL
        assert new_domain.created_by_antismash
        assert new_domain.locus_tag == "locus"
        assert new_domain.protein_location == self.protein_location

    def test_subtype_recognition(self):
        assert self.domain.domain_id
        bio = self.domain.to_biopython()[0]
        assert bio.qualifiers["aSTool"] == [TOOL]
        assert bio.qualifiers["domain_id"]
        new = AntismashDomain.from_biopython(bio)
        assert isinstance(new, ModularDomain)

    def test_linebreaks(self):
        self.domain.domain_id = "some-long-name.Condensation.1"
        bio = self.domain.to_biopython()
        # pretend a genbank was written and is being read back in
        # but the line was long enough to cause a line break
        bio[0].qualifiers["domain_id"][0] = self.domain.domain_id.replace(".", ". ")
        new = ModularDomain.from_biopython(bio[0])
        assert new.domain_id == self.domain.domain_id
