# License: GNU Affero General Public License v3 or later
# A copy of GNU AGPL v3 should have been included in this software package in LICENSE.txt.

""" Classes representing qualifiers for type II PKS details.
"""

from typing import Dict, List, Optional

from .secmet import _parse_format


class T2PKSQualifier:
    """A qualifier for type II PKS details at a Cluster level.
    """
    STARTER_KEY = "t2pks_starter_units"
    MAL_ELONGATION_KEY = "t2pks_malonyl_elongations"
    WEIGHTS_KEY = "t2pks_molecular_weights"
    PRODUCT_CLASSES_KEY = "t2pks_product_classes"

    WEIGHT_TEMPLATE = "{} (Da): {:.3f}"

    def __init__(self, starter_units: List[str], malonyl_elongations: List[str],
                 product_classes: List[str], molecular_weights: Dict[str, float]) -> None:
        if not starter_units:
            raise ValueError("cannot have a type II PKS qualifier with no starter unit")
        if bool(malonyl_elongations) != bool(molecular_weights):
            raise ValueError("both malonyl elongations and molecular weights must be supplied if one is")

        self.starter_units = starter_units
        self.malonyl_elongations = malonyl_elongations
        self.product_classes = product_classes
        self.molecular_weights = molecular_weights

    def to_biopython_qualifiers(self) -> Dict[str, List[str]]:
        """ Converts the qualifier into a BioPython-compatible format """
        qualifiers = {self.STARTER_KEY: self.starter_units}
        if self.malonyl_elongations:
            qualifiers[self.MAL_ELONGATION_KEY] = self.malonyl_elongations
            weights = []
            for combination, weight in self.molecular_weights.items():
                weights.append(self.WEIGHT_TEMPLATE.format(combination, weight))
            qualifiers[self.WEIGHTS_KEY] = weights
        if self.product_classes:
            qualifiers[self.PRODUCT_CLASSES_KEY] = self.product_classes
        return qualifiers

    @staticmethod
    def from_biopython_qualifiers(qualifiers: Dict[str, List[str]]) -> Optional["T2PKSQualifier"]:
        """ Converts the relevant section of a BioPython qualifiers dictionary
            to a T2PKSQualifier.

            Removes the relevant fields from the input dictionary.

            Arguments:
                qualifiers: a list of strings

            Returns:
                A T2PKSQualifier instance with the same details as the provided
                qualifier.

        """
        starters = qualifiers.pop(T2PKSQualifier.STARTER_KEY, [])
        if not starters:
            return None

        elongations = qualifiers.pop(T2PKSQualifier.MAL_ELONGATION_KEY, [])
        raw_weights = qualifiers.pop(T2PKSQualifier.WEIGHTS_KEY, [])
        classes = qualifiers.pop(T2PKSQualifier.PRODUCT_CLASSES_KEY, [])

        weights = {}
        for raw_weight in raw_weights:
            combination, weight = _parse_format(T2PKSQualifier.WEIGHT_TEMPLATE, raw_weight)
            weights[combination] = float(weight)
        return T2PKSQualifier(starters, elongations, classes, weights)
