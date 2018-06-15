# License: GNU Affero General Public License v3 or later
# A copy of GNU AGPL v3 should have been included in this software package in LICENSE.txt.

""" Classes representing complex qualifiers for features.
"""

from typing import Dict, List


class GOQualifier:
    """A qualifier for tracking Gene Ontology terms for a PFAM domain.
        Cannot be directly used as a qualifier for BioPython's SeqFeature.
    """
    def __init__(self, go_entries: Dict[str, str]) -> None:  # dict mapping Gene Ontology IDs to readable descriptions
        self.go_entries = go_entries
        self.ids = list(go_entries)
        self.descriptions = list(go_entries.values())

    def to_biopython(self) -> List[str]:
        """Convert GOQualifier to BioPython-style qualifier."""
        return ["{}: {}".format(go_id, go_description) for go_id, go_description in sorted(self.go_entries.items())]

    @staticmethod
    def from_biopython(qualifier: List[str]) -> "GOQualifier":
        """Convert BioPython-style qualifier to GOQualifier.

            Arguments:
                qualifier: BioPython-style qualifier (list of strings)

            Returns:
                A GOQualifier instance constructed from the qualifier.

        """
        go_entries = {}
        for go_string in qualifier:
            go_id, separator, go_description = go_string.partition(": ")
            if not separator:
                raise ValueError("Cannot parse qualifier: %s" % qualifier)
            go_entries[go_id] = go_description
        result = GOQualifier(go_entries)
        return result
