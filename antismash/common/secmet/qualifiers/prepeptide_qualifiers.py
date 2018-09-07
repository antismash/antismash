# License: GNU Affero General Public License v3 or later
# A copy of GNU AGPL v3 should have been included in this software package in LICENSE.txt.

""" Qualifiers to contain extra annotations for RiPP prepeptides
"""

from typing import Dict, List, Optional
from typing import Type  # comment hints, pylint: disable=unused-import


class RiPPQualifier:
    """ A generic qualifier for RiPP annotations """
    __slots__ = ["rodeo_score"]

    def __init__(self, rodeo_score: int = 0) -> None:
        self.rodeo_score = rodeo_score

    def to_biopython_qualifiers(self) -> Dict[str, List[str]]:
        """ Generate a biopython-like dictionary of qualifiers containing the
            information from this qualifier
        """
        qualifiers = {}
        if self.rodeo_score is not 0:
            qualifiers["RODEO_score"] = [str(self.rodeo_score)]
        return qualifiers

    @classmethod
    def from_biopython_qualifiers(cls, qualifiers: Dict[str, List[str]]) -> "RiPPQualifier":
        """ Rebuilds an instance of this qualifier from a biopython-like dictionary
            of qualifiers. Removes the relevant sections of the qualifiers. """
        return cls(int(qualifiers.pop("RODEO_score")[0]))


class LanthiQualifier(RiPPQualifier):
    """ A qualifier for lanthipeptide-specific annotations """
    __slots__ = ["lan_bridges", "aminovinyl_group", "chlorinated", "oxygenated", "lactonated"]

    def __init__(self, lan_bridges: int, rodeo_score: int,  # pylint: disable=too-many-arguments
                 aminovinyl_group: bool, chlorinated: bool, oxygenated: bool,
                 lactonated: bool) -> None:
        super().__init__(rodeo_score)
        self.lan_bridges = lan_bridges
        self.aminovinyl_group = aminovinyl_group
        self.chlorinated = chlorinated
        self.oxygenated = oxygenated
        self.lactonated = lactonated

    def get_modifications(self) -> List[str]:
        """ Returns the various modifications of the lanthipeptide
        """
        mods = []
        if self.aminovinyl_group:
            mods.append("AviCys")
        if self.chlorinated:
            mods.append("Cl")
        if self.oxygenated:
            mods.append("OH")
        if self.lactonated:
            mods.append("Lac")
        return mods

    def to_biopython_qualifiers(self) -> Dict[str, List[str]]:
        qualifiers = super().to_biopython_qualifiers()
        qualifiers["number_of_bridges"] = [str(self.lan_bridges)]
        mods = self.get_modifications()
        if mods:
            qualifiers["predicted_additional_modifications"] = mods
        return qualifiers

    @classmethod
    def from_biopython_qualifiers(cls, qualifiers: Dict[str, List[str]]) -> "LanthiQualifier":
        mods = qualifiers.pop("predicted_additional_modifications", [])
        aminovinyl = "AviCys" in mods
        chlorinated = "Cl" in mods
        oxygenated = "OH" in mods
        lactonated = "Lac" in mods
        return cls(int(qualifiers.pop("number_of_bridges")[0]),
                   int(qualifiers.pop("RODEO_score")[0]),
                   aminovinyl, chlorinated, oxygenated, lactonated)

def rebuild_qualifier(data: Dict[str, List[str]], kind: str) -> Optional[RiPPQualifier]:
    """ Rebuilds a relevant RiPPQualifier for the given kind from the provided
        biopython qualifiers.

        Removes used portions the qualifiers.

        Arguments:
            data: the qualifiers in biopython format
            kind: the peptide class

        Returns:
            a RiPPQualifier subclass matching the peptide class provided
    """
    if not data or "RODEO_score" not in data:
        return None
    classes = {
        "lanthipeptide": LanthiQualifier,
    }  # type: Dict[str, Type[RiPPQualifier]]
    if kind not in classes:
        raise ValueError("no known qualifier builder for prepeptide kind: %s" % kind)
    return classes[kind].from_biopython_qualifiers(data)
