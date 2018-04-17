# License: GNU Affero General Public License v3 or later
# A copy of GNU AGPL v3 should have been included in this software package in LICENSE.txt.

# for test files, silence irrelevant and noisy pylint warnings
# pylint: disable=no-self-use,protected-access,missing-docstring

import unittest

from antismash.common.hmmscan_refinement import HMMResult
from antismash.common.secmet.qualifiers import NRPSPKSQualifier, GOQualifier

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
        for go_id in go_qualifier.ids:
            assert go_id in original_go_entries

    def test_go_descs(self):
        original_go_entries = {'GO:0004871': 'signal transducer activity', 'GO:0007165': 'signal transduction',
                               'GO:0016020': 'membrane'}
        go_qualifier = GOQualifier(original_go_entries)
        for go_description in go_qualifier.descriptions:
            assert go_description in original_go_entries.values()

    def test_biopython_to_and_from(self):
        original = GOQualifier({'GO:0004871': 'signal transducer activity', 'GO:0007165': 'signal transduction',
                                'GO:0016020': 'membrane'})
        new = GOQualifier.from_biopython(original.to_biopython())
        assert original.go_entries == new.go_entries

    def test_parse_broken_qualifier(self):
        broken_qualifier = ["GO:0004871: signal transducer activity", "GO:0007165; signal transduction"]
        with self.assertRaisesRegex(ValueError, "Cannot parse qualifier"):
            GOQualifier.from_biopython(broken_qualifier)
