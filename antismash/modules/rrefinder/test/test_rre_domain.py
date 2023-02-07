# License: GNU Affero General Public License v3 or later
# A copy of GNU AGPL v3 should have been included in this software package in LICENSE.txt.

# for test files, silence irrelevant and noisy pylint warnings
# pylint: disable=use-implicit-booleaness-not-comparison,protected-access,missing-docstring

import unittest

from antismash.common.secmet.features import AntismashDomain, FeatureLocation

from antismash.modules.rrefinder.rre_domain import RREDomain


class TestRRE(unittest.TestCase):
    def setUp(self):
        self.protein_location = FeatureLocation(1, 5)
        self.location = FeatureLocation(6, 10)
        self.domain = 'RRE_type_a'
        self.description = 'This is a test RRE'
        self.locus_tag = 'locus_tag_a'
        self.identifier = 'RREFam001'
        self.version = 1
        self.full_identifier = f"{self.identifier}.{self.version}"
        self.rre = RREDomain(self.location, self.description, self.protein_location,
                             self.full_identifier, self.locus_tag, self.domain)
        self.rre.domain_id = f"{self.locus_tag}_{self.identifier}_1"

    def test_init(self):
        assert self.rre.locus_tag == self.locus_tag
        assert self.rre.description == self.description
        assert self.rre.domain == self.domain
        assert self.rre.location == self.location
        assert self.rre.protein_location == self.protein_location
        assert self.rre.identifier == self.identifier
        assert self.rre.version == self.version

    def test_init_no_version(self):
        with self.assertRaisesRegex(ValueError, "missing version"):
            RREDomain(self.location, self.description, self.protein_location,
                      self.identifier, self.locus_tag, self.domain)

    def test_init_wrong_description(self):
        # Test wrong description entries
        with self.assertRaises(TypeError):
            RREDomain(self.location, 5, self.protein_location, self.identifier,
                      self.locus_tag, self.domain)
        with self.assertRaisesRegex(ValueError, "RRE description cannot be empty"):
            RREDomain(self.location, '', self.protein_location, self.identifier,
                      self.locus_tag, self.domain)
        with self.assertRaisesRegex(ValueError, "RREFam identifier cannot be empty"):
            RREDomain(self.location, self.description, self.protein_location, '',
                      self.locus_tag, self.domain)
        with self.assertRaises(ValueError):
            RREDomain(self.location, self.description, self.protein_location, 'not_a_valid_identifier',
                      self.locus_tag, self.domain)

    def test_conversion(self):
        bio = self.rre.to_biopython()

        assert len(bio) == 1
        assert bio[0].qualifiers['description'][0] == self.description
        assert bio[0].qualifiers['identifier'][0] == self.full_identifier

        rre = RREDomain.from_biopython(bio[0])
        assert rre.description == self.rre.description
        assert rre.protein_location == self.rre.protein_location
        assert rre.locus_tag == self.rre.locus_tag
        assert rre.identifier == self.rre.identifier
        assert rre.version == self.rre.version

        # Test with extra qualifiers
        bio = self.rre.to_biopython(qualifiers={'some_qualifier': ['some_value']})
        assert bio[0].qualifiers.get('some_qualifier')
        assert bio[0].qualifiers['some_qualifier'][0] == 'some_value'

    def test_subvariant_detection(self):
        bio = self.rre.to_biopython()[0]
        new = AntismashDomain.from_biopython(bio)
        assert isinstance(new, RREDomain)
