# License: GNU Affero General Public License v3 or later
# A copy of GNU AGPL v3 should have been included in this software package in LICENSE.txt.

""" A qualifier for ActiveSiteFinder annotations """

from typing import List
from typing import Set  # comment hint  # pylint: disable=unused-import


class ActiveSiteFinderQualifier:
    """ A qualifier for tracking active sites found (or not found) within a CDS.
    """
    def __init__(self) -> None:
        self._hits = set()  # type: Set[str]

    @property
    def hits(self) -> List[str]:
        """ Returns a list of all active site notes """
        return sorted(list(self._hits))

    def add(self, label: str) -> None:
        """ Adds an active site presence label to the qualifier. """
        self._hits.add(str(label))

    def to_biopython(self) -> List[str]:
        """ Creates a BioPython style qualifier from the qualifier """
        return self.hits

    def __bool__(self) -> bool:
        return bool(self._hits)
