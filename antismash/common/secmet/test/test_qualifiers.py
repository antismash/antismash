# License: GNU Affero General Public License v3 or later
# A copy of GNU AGPL v3 should have been included in this software package in LICENSE.txt.

# for test files, silence irrelevant and noisy pylint warnings
# pylint: disable=no-self-use,protected-access,missing-docstring

import unittest

from antismash.common.hmmscan_refinement import HMMResult
from antismash.common.secmet.feature import CDSFeature, FeatureLocation
from antismash.common.secmet.qualifiers import NRPSPKSQualifier, GOQualifier, GeneFunction


class TestNRPSPKS(unittest.TestCase):
    def test_counter(self):
        qualifier = NRPSPKSQualifier(strand=1)
        types = [("PKS_AT", "_AT"), ("PKS_KR", "_KR"), ("CAL_domain", "_CAL"),
                 ("AMP-binding", "_A"), ("PKS_KS", "_KS"), ("ACP", "_OTHER")]
        expected = set()
        for pks_type, suffix in types:
            domain = HMMResult(pks_type, 1, 1, 1, 1)
            suffix = suffix + "%d"
            for i in range(3):
                qualifier.add_domain(domain, "missing")
                expected.add(suffix % (i + 1))
        assert len(qualifier.domains) == 3 * len(types)
        assert {domain.label for domain in qualifier.domains} == expected

    def test_no_append(self):
        qualifier = NRPSPKSQualifier(strand=1)
        with self.assertRaisesRegex(NotImplementedError, "Appending to this list won't work"):
            qualifier.append("test")

        with self.assertRaisesRegex(NotImplementedError, "Extending this list won't work"):
            qualifier.extend(["test"])

    def test_biopython_compatibility(self):
        qualifier = NRPSPKSQualifier(strand=1)
        assert isinstance(qualifier, list)
        for pks in ["PKS_AT", "AMP-binding"]:
            qualifier.add_domain(HMMResult(pks, 1, 1, 1, 1), "missing")
            qualifier.add_subtype(pks + "dummy")
        assert len(qualifier) == 4
        for i in qualifier:
            assert isinstance(i, str)


class TestGeneFunction(unittest.TestCase):
    def test_membership(self):
        assert GeneFunction.OTHER
        with self.assertRaises(AttributeError):
            print(GeneFunction.non_existant)

    def test_equality(self):
        assert GeneFunction.OTHER == GeneFunction.OTHER
        assert GeneFunction.CORE != GeneFunction.OTHER

    def test_string_conversion(self):
        assert str(GeneFunction.CORE) == "biosynthetic"
        assert str(GeneFunction.ADDITIONAL) == "biosynthetic-additional"
        assert str(GeneFunction.OTHER) == "other"
        for member in dir(GeneFunction):
            if not member.isupper() or member in ["CORE", "ADDITIONAL"]:
                continue
            assert str(getattr(GeneFunction, member)) == member.lower()

    def test_cds_function(self):
        cds = CDSFeature(FeatureLocation(1, 5, 1), locus_tag="foo")
        # default value
        assert cds.gene_functions.get_classification() == GeneFunction.OTHER
        assert cds.gene_function == GeneFunction.OTHER
        # check bad values can't be assigned
        with self.assertRaises(AssertionError):
            cds.gene_functions.add("other", "a", "b")
        with self.assertRaises(AttributeError):
            cds.gene_functions = 0

        cds.gene_functions.add(GeneFunction.ADDITIONAL, "first_tool", "dummy")
        assert cds.gene_functions.get_classification() == GeneFunction.ADDITIONAL
        assert cds.gene_function == GeneFunction.ADDITIONAL
        # conflicting, so back to OTHER
        cds.gene_functions.add(GeneFunction.TRANSPORT, "other_tool", "dummy")
        assert cds.gene_functions.get_classification() == GeneFunction.OTHER
        assert cds.gene_function == GeneFunction.OTHER
        # but smcogs overrides that
        cds.gene_functions.add(GeneFunction.REGULATORY, "smcogs", "dummy")
        assert cds.gene_functions.get_classification() == GeneFunction.REGULATORY
        # and cluster definition overrides even that
        cds.gene_functions.add(GeneFunction.CORE, "cluster_definition", "dummy")
        assert cds.gene_functions.get_classification() == GeneFunction.CORE

        # and that we still have tracked these
        smcogs = cds.gene_functions.get_by_tool("smcogs")
        assert len(smcogs) == 1
        assert smcogs[0].function == GeneFunction.REGULATORY

        adds = cds.gene_functions.get_by_function(GeneFunction.ADDITIONAL)
        assert len(adds) == 1
        assert adds[0].tool == "first_tool"

    def test_cds_function_conversion(self):
        cds = CDSFeature(FeatureLocation(1, 5, 1), locus_tag="foo")
        assert cds.gene_function == GeneFunction.OTHER
        assert CDSFeature.from_biopython(cds.to_biopython()[0]).gene_function == GeneFunction.OTHER
        cds.gene_functions.add(GeneFunction.ADDITIONAL, "tool", "desc")
        assert cds.gene_function == GeneFunction.ADDITIONAL
        assert CDSFeature.from_biopython(cds.to_biopython()[0]).gene_function == GeneFunction.ADDITIONAL


class TestGOQualifier(unittest.TestCase):
    def test_go_entries(self):
        original_go_entries = {'GO:0004871': 'signal transducer activity', 'GO:0007165': 'signal transduction',
                               'GO:0016020': 'membrane'}
        go_qualifier = GOQualifier(original_go_entries)
        assert go_qualifier.go_entries == original_go_entries

    def test_go_ids(self):
        original_go_entries = {'GO:0004871': 'signal transducer activity', 'GO:0007165': 'signal transduction',
                               'GO:0016020': 'membrane'}
        go_qualifier = GOQualifier(original_go_entries)
        assert set(go_qualifier.ids) == set(original_go_entries)

    def test_go_descs(self):
        original_go_entries = {'GO:0004871': 'signal transducer activity', 'GO:0007165': 'signal transduction',
                               'GO:0016020': 'membrane'}
        go_qualifier = GOQualifier(original_go_entries)
        assert set(go_qualifier.descriptions) == set(original_go_entries.values())

    def test_biopython_to_and_from(self):
        original = GOQualifier({'GO:0004871': 'signal transducer activity', 'GO:0007165': 'signal transduction',
                                'GO:0016020': 'membrane'})
        new = GOQualifier.from_biopython(original.to_biopython())
        assert original.go_entries == new.go_entries

    def test_parse_broken_qualifier(self):
        # test if wrong separator between GO ID and description (semicolon instead of colon) is caught
        broken_qualifier = ["GO:0004871: signal transducer activity", "GO:0007165; signal transduction"]
        with self.assertRaisesRegex(ValueError, "Cannot parse qualifier"):
            GOQualifier.from_biopython(broken_qualifier)
