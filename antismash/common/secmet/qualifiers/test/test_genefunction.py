# License: GNU Affero General Public License v3 or later
# A copy of GNU AGPL v3 should have been included in this software package in LICENSE.txt.

# for test files, silence irrelevant and noisy pylint warnings
# pylint: disable=no-self-use,protected-access,missing-docstring

import unittest

from ...features import CDSFeature, FeatureLocation
from ..gene_functions import GeneFunction


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
        cds = CDSFeature(FeatureLocation(1, 5, 1), locus_tag="foo", translation="A")
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
        cds.gene_functions.add(GeneFunction.CORE, "cluster_definition", "dummy", "product")
        assert cds.gene_functions.get_classification() == GeneFunction.CORE

        # and that we still have tracked these
        smcogs = cds.gene_functions.get_by_tool("smcogs")
        assert len(smcogs) == 1
        assert smcogs[0].function == GeneFunction.REGULATORY

        adds = cds.gene_functions.get_by_function(GeneFunction.ADDITIONAL)
        assert len(adds) == 1
        assert adds[0].tool == "first_tool"

    def test_cds_function_conversion(self):
        cds = CDSFeature(FeatureLocation(1, 5, 1), locus_tag="foo", translation="A")
        assert cds.gene_function == GeneFunction.OTHER
        assert CDSFeature.from_biopython(cds.to_biopython()[0]).gene_function == GeneFunction.OTHER
        cds.gene_functions.add(GeneFunction.ADDITIONAL, "tool", "desc")
        assert cds.gene_function == GeneFunction.ADDITIONAL
        assert CDSFeature.from_biopython(cds.to_biopython()[0]).gene_function == GeneFunction.ADDITIONAL
