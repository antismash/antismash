# License: GNU Affero General Public License v3 or later
# A copy of GNU AGPL v3 should have been included in this software package in LICENSE.txt.

""" A class for CDS motif features """

from collections import OrderedDict
from typing import Dict, List, Optional

from Bio.SeqFeature import SeqFeature

from .domain import Domain
from .feature import Feature, Location


class CDSMotif(Domain):
    """ A base class for features that represent a motif within a CDSFeature """
    __slots__ = ["motif"]

    def __init__(self, location: Location, tool: Optional[str] = None) -> None:
        # if there's a tool, it was created by antismash
        created = tool is not None
        super().__init__(location, feature_type="CDS_motif", tool=tool, created_by_antismash=created)

    @staticmethod
    def from_biopython(bio_feature: SeqFeature, feature: Optional["CDSMotif"] = None,  # type: ignore
                       leftovers: Optional[Dict[str, List[str]]] = None) -> "CDSMotif":
        if leftovers is None:
            leftovers = Feature.make_qualifiers_copy(bio_feature)
        if not feature:
            feature = CDSMotif(bio_feature.location, leftovers.pop("aSTool", [""])[0] or None)

        updated = super(CDSMotif, feature).from_biopython(bio_feature, feature, leftovers)
        assert updated is feature
        assert isinstance(updated, CDSMotif)
        return updated

    def to_biopython(self, qualifiers: Dict[str, List] = None) -> List[SeqFeature]:
        mine = OrderedDict()  # type: Dict[str, List[str]]
        if qualifiers:
            mine.update(qualifiers)
        return super().to_biopython(mine)
