# License: GNU Affero General Public License v3 or later
# A copy of GNU AGPL v3 should have been included in this software package in LICENSE.txt.

""" A class for subregion features """

from typing import Dict, List, Optional

from Bio.SeqFeature import SeqFeature

from .cdscollection import CDSCollection
from .feature import FeatureLocation, Feature


class SubRegion(CDSCollection):
    """ A feature which marks a specific region of a record as interesting,
        without being considered a cluster.
    """
    def __init__(self, location: FeatureLocation, tool: str, probability: float = None, anchor: str = "") -> None:
        super().__init__(location, feature_type="subregion")
        self.tool = tool
        self.probability = probability
        self.anchor = anchor  # if anchored to a gene/CDS, this is the name

    def get_subregion_number(self) -> int:
        """ Returns the subregion's numeric ID, only guaranteed to be consistent
            when the same subregions are defined in the parent record
        """
        if not self._parent_record:
            raise ValueError("SubRegion not in a record")
        return self._parent_record.get_subregion_number(self)

    def to_biopython(self, qualifiers: Optional[Dict[str, List[str]]] = None) -> List[SeqFeature]:
        if qualifiers is None:
            qualifiers = {}
        if self._parent_record:
            qualifiers["subregion_number"] = [str(self.get_subregion_number())]
        qualifiers["aStool"] = [self.tool]
        if self.probability is not None:
            qualifiers["probability"] = [str(self.probability)]
        if self.anchor:
            qualifiers["anchor"] = [self.anchor]
        return super().to_biopython(qualifiers)

    @staticmethod
    def from_biopython(bio_feature: SeqFeature, feature: "SubRegion" = None,  # type: ignore
                       leftovers: Optional[Dict] = None) -> "SubRegion":
        if leftovers is None:
            leftovers = Feature.make_qualifiers_copy(bio_feature)

        tool = leftovers.pop("aStool")[0]
        probability = None
        if "probability" in leftovers:
            probability = float(leftovers.pop("probability")[0])
        anchor = leftovers.pop("anchor", [""])[0]
        if not feature:
            feature = SubRegion(bio_feature.location, tool, probability, anchor)

        # remove the subregion_number, as it's not relevant
        leftovers.pop("subregion_number", "")

        # grab parent optional qualifiers
        super(SubRegion, feature).from_biopython(bio_feature, feature=feature, leftovers=leftovers)
        return feature
