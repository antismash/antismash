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
        if not isinstance(rodeo_score, int):
            raise TypeError("RODEO score must be an int, not %s" % type(rodeo_score))
        self.rodeo_score = rodeo_score

    def to_biopython_qualifiers(self) -> Dict[str, List[str]]:
        """ Generate a biopython-like dictionary of qualifiers containing the
            information from this qualifier
        """
        return {"RODEO_score": [str(self.rodeo_score)]}

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


class ThioQualifier(RiPPQualifier):
    """ A qualifier for thiopeptide-specific annotations """
    __slots__ = ["amidation", "macrocycle", "core_features", "mature_weights"]

    def __init__(self, rodeo_score: int, amidation: bool, macrocycle: str,  # pylint: disable=too-many-arguments
                 core_features: str, mature_weights: List[float]) -> None:
        super().__init__(rodeo_score)
        self.amidation = amidation
        self.macrocycle = macrocycle
        self.core_features = core_features
        self.mature_weights = mature_weights

    @property
    def tail_reaction(self) -> str:
        """ The tail reaction of a thiopeptide if present """
        if self.amidation:
            return "dealkylation of C-Terminal residue; amidation"
        return ''

    def to_biopython_qualifiers(self) -> Dict[str, List[str]]:
        qualifiers = super().to_biopython_qualifiers()
        qualifiers.update({
            "macrocycle": [self.macrocycle],
            "core_features": [self.core_features],
            "mature_weights": [str(weight) for weight in self.mature_weights],
        })
        if self.amidation:
            qualifiers["tail_reaction"] = [self.tail_reaction]
        return qualifiers

    @classmethod
    def from_biopython_qualifiers(cls, qualifiers: Dict[str, List[str]]) -> "ThioQualifier":
        weights = [float(weight) for weight in qualifiers.pop("mature_weights")]
        return cls(int(qualifiers.pop("RODEO_score")[0]),
                   qualifiers.pop("tail_reaction", [""])[0].endswith("amidation"),
                   qualifiers.pop("macrocycle")[0],
                   qualifiers.pop("core_features")[0],
                   weights)


class LassoQualifier(RiPPQualifier):
    """ A qualifier for lassopeptide-specific annotations """
    def __init__(self, rodeo_score: int, num_bridges: int,  # pylint: disable=too-many-arguments
                 macrolactam: str, cut_mass: float, cut_weight: float) -> None:
        super().__init__(rodeo_score)
        self.num_bridges = num_bridges
        self.macrolactam = macrolactam
        self.cut_mass = cut_mass
        self.cut_weight = cut_weight

    def to_biopython_qualifiers(self) -> Dict[str, List[str]]:
        qualifiers = super().to_biopython_qualifiers()
        qualifiers.update({
            "number_of_bridges": [str(self.num_bridges)],
            "macrolactam": [self.macrolactam],
            "cut_mass": [str(self.cut_mass)],
            "cut_weight": [str(self.cut_weight)],
        })
        return qualifiers

    @classmethod
    def from_biopython_qualifiers(cls, qualifiers: Dict[str, List[str]]) -> "LassoQualifier":
        return cls(int(qualifiers.pop("RODEO_score")[0]),
                   int(qualifiers.pop("number_of_bridges")[0]),
                   qualifiers.pop("macrolactam")[0],
                   float(qualifiers.pop("cut_mass")[0]),
                   float(qualifiers.pop("cut_weight")[0]))


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
        "thiopeptide": ThioQualifier,
        "lassopeptide": LassoQualifier,
    }  # type: Dict[str, Type[RiPPQualifier]]
    if kind not in classes:
        raise ValueError("no known qualifier builder for prepeptide kind: %s" % kind)
    return classes[kind].from_biopython_qualifiers(data)
