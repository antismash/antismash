# License: GNU Affero General Public License v3 or later
# A copy of GNU AGPL v3 should have been included in this software package in LICENSE.txt.

# for test files, silence irrelevant and noisy pylint warnings
# pylint: disable=no-self-use,protected-access,missing-docstring

import unittest

from antismash.common.secmet.qualifiers import GOQualifier


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
