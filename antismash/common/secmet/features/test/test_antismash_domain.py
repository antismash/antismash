# License: GNU Affero General Public License v3 or later
# A copy of GNU AGPL v3 should have been included in this software package in LICENSE.txt.

# for test files, silence irrelevant and noisy pylint warnings
# pylint: disable=no-self-use,protected-access,missing-docstring

import unittest

from antismash.common.secmet.features import AntismashDomain, Feature, FeatureLocation
from antismash.common.secmet.features import antismash_domain


class TestConversion(unittest.TestCase):
    def setUp(self):
        self.protein_location = FeatureLocation(0, 1)
        self.domain = AntismashDomain(FeatureLocation(1, 3, 1), locus_tag="locus",
                                      tool="test", protein_location=self.protein_location)
        self.domain.domain_id = "test1"

    def test_conversion(self):
        domain = self.domain
        domain.asf.add("first")
        domain.asf.add("second")
        assert domain.tool == "test"
        assert domain.created_by_antismash
        assert domain.locus_tag == "locus"

        bio = domain.to_biopython()
        assert len(bio) == 1
        assert bio[0].qualifiers["aSTool"] == ["test"]
        assert bio[0].qualifiers["tool"] == ["antismash"]
        new_domain = AntismashDomain.from_biopython(bio[0])
        assert new_domain.asf.hits == domain.asf.hits
        assert new_domain.asf.hits == ["first", "second"]
        assert new_domain.tool == domain.tool == "test"
        assert new_domain.created_by_antismash
        assert new_domain.locus_tag == "locus"
        assert new_domain.protein_location == self.protein_location

    def test_linebreaks(self):
        self.domain.domain_id = "some-long-name.Condensation.1"
        bio = self.domain.to_biopython()
        # pretend a genbank was written and is being read back in
        # but the line was long enough to cause a line break
        bio[0].qualifiers["domain_id"][0] = self.domain.domain_id.replace(".", ". ")
        new = AntismashDomain.from_biopython(bio[0])
        assert new.domain_id == self.domain.domain_id


class TestSubtyping(unittest.TestCase):
    def setUp(self):
        self.original = antismash_domain._SUBTYPE_MAPPING
        antismash_domain._SUBTYPE_MAPPING = {}

    def tearDown(self):
        antismash_domain._SUBTYPE_MAPPING = self.original

    def register(self, name, subclass):
        antismash_domain.register_asdomain_variant(name, subclass)

    def test_invalid_class(self):
        for value in [int, Feature]:
            with self.assertRaisesRegex(TypeError, "not a subclass"):
                self.register("test", value)

    def test_duplicate_tools(self):
        self.register("test", AntismashDomain)
        with self.assertRaisesRegex(ValueError, "already present as a subtype"):
            self.register("test", AntismashDomain)


class TestValues(unittest.TestCase):
    def test_domain_passthrough(self):
        for value in ["test", "DomainName"]:
            domain = AntismashDomain(FeatureLocation(1, 3, 1), locus_tag="locus",
                                     tool="test", protein_location=FeatureLocation(0, 1),
                                     domain=value)
            assert domain.domain == value
