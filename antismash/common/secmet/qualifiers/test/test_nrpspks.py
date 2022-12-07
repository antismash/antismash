# License: GNU Affero General Public License v3 or later
# A copy of GNU AGPL v3 should have been included in this software package in LICENSE.txt.

# for test files, silence irrelevant and noisy pylint warnings
# pylint: disable=use-implicit-booleaness-not-comparison,protected-access,missing-docstring

import unittest

from antismash.common.hmmscan_refinement import HMMResult

from ..nrps_pks import NRPSPKSQualifier


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

    def test_biopython_compatibility(self):
        qualifier = NRPSPKSQualifier(strand=1)
        for pks in ["PKS_AT", "AMP-binding"]:
            qualifier.add_domain(HMMResult(pks, 1, 1, 1, 1), "missing")
        assert len(qualifier) == 2
        for i in qualifier:
            assert isinstance(i, str)

    def test_biopython_conversion(self):
        qualifier = NRPSPKSQualifier(strand=1)
        for name in ["PKS_AT", "AMP-binding"]:
            domain = qualifier.add_domain(HMMResult(name, 1, 1, 1, 1), "missing")
            domain.subtypes = [name + "sub"]
        qualifier.type = "some type"

        bio = list(qualifier)
        for val in bio:
            assert isinstance(val, str)

        new = NRPSPKSQualifier(strand=1)
        new.add_from_qualifier(bio)
        for old_domain, new_domain in zip(qualifier.domains, new.domains):
            assert old_domain.name == new_domain.name
            assert old_domain.subtypes == new_domain.subtypes
        assert list(qualifier) == list(new)

        for bad in [["mismatching info"], ["Domain: missing info"]]:
            with self.assertRaisesRegex(ValueError, "unknown NRPS/PKS qualifier|could not match"):
                new.add_from_qualifier(bad)
