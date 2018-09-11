# License: GNU Affero General Public License v3 or later
# A copy of GNU AGPL v3 should have been included in this software package in LICENSE.txt.

# for test files, silence irrelevant and noisy pylint warnings
# pylint: disable=no-self-use,protected-access,missing-docstring

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
