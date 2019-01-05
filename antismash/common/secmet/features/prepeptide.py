# License: GNU Affero General Public License v3 or later
# A copy of GNU AGPL v3 should have been included in this software package in LICENSE.txt.

""" A class for representing RiPP prepeptides """

from typing import Any, Dict, List
from typing import Optional  # comment hints, pylint: disable=unused-import

from Bio.SeqFeature import SeqFeature

from ..errors import SecmetInvalidInputError
from .cds_motif import CDSMotif
from .feature import Feature, Location
from ..locations import build_location_from_others, location_from_string
from ..qualifiers.prepeptide_qualifiers import RiPPQualifier  # comment hints, pylint: disable=unused-import
from ..qualifiers.prepeptide_qualifiers import rebuild_qualifier


class Prepeptide(CDSMotif):  # pylint: disable=too-many-instance-attributes
    """ A class representing a prepeptide. Used for tracking a multi-feature
        construction with a leader, core and tail. To allow for multiple types
        of prepeptide (e.g. lanthi- or sacti-peptides), only the core must exist.
    """
    def __init__(self, location: Location, peptide_class: str, core: str, locus_tag: str,
                 tool: str, peptide_subclass: str = None, score: float = 0., monoisotopic_mass: float = 0.,
                 molecular_weight: float = 0., alternative_weights: List[float] = None,
                 leader: str = "", tail: str = "") -> None:
        """
            Arguments:
                peptide_class: the kind of prepeptide, e.g. 'lanthipeptide', 'thiopeptide'
                core: the sequence of the core
                locus_tag: the locus tag to use for the feature
                tool: the name of the tool responsible for creating the prepeptide
                prepeptide_subclass: the subclass of the prepeptide, e.g. 'Type II'
                score: the prepeptide score
                monoisotopic_mass: the monoisotopic mass of the prepeptide
                molecular_weight: the molecular weight of the prepeptide
                alternative_weights: a list of alternative weights for the prepeptide
                leader: the sequence of the leader, if it exists
                tail: the sequence of the tail, if it exists
        """
        if not tool:
            raise ValueError("tool must not be empty or None")
        super().__init__(location, tool=tool)
        for arg in [peptide_class, core, leader, tail]:
            assert isinstance(arg, str), type(arg)
        self._leader = leader
        if not core:
            raise ValueError("Prepeptides must have a core")
        self._core = core
        self._tail = tail
        self.locus_tag = locus_tag
        self.peptide_class = peptide_class
        if peptide_subclass:
            peptide_subclass = peptide_subclass.replace("-", " ")  # "Type-II" > "Type II"
        self.peptide_subclass = peptide_subclass
        self.score = float(score)
        self.monoisotopic_mass = float(monoisotopic_mass)
        self.molecular_weight = float(molecular_weight)
        self.alternative_weights = []  # type: List[float]
        if alternative_weights is not None:
            self.alternative_weights = [float(weight) for weight in alternative_weights]

        self.detailed_information = None  # type: Optional[RiPPQualifier]

    @property
    def translation(self) -> str:
        return self._leader + self._core + self._tail

    @translation.setter
    def translation(self) -> None:
        raise AttributeError("Cannot assign to translation in a Prepeptide")

    @property
    def leader(self) -> str:
        """ The leader sequence of the prepeptide """
        return self._leader

    @leader.setter
    def leader(self, leader: str) -> None:
        assert isinstance(leader, str)
        self._leader = leader

    @property
    def core(self) -> str:
        """ The core sequence of the prepeptide """
        return self._core

    @core.setter
    def core(self, core: str) -> None:
        assert isinstance(core, str)
        self._core = core

    @property
    def tail(self) -> str:
        """ The tail sequence of the prepeptide """
        return self._tail

    @tail.setter
    def tail(self, tail: str) -> None:
        assert isinstance(tail, str)
        self._tail = tail

    def get_name(self) -> str:
        """ Returns the locus tag of the parent CDS.

            Uses the same function name as the CDSFeature for consistency.
        """
        assert isinstance(self.locus_tag, str) and self.locus_tag
        return self.locus_tag

    def to_biopython(self, qualifiers: Dict[str, List] = None) -> List[SeqFeature]:
        """ Generates up to three SeqFeatures, depending if leader and tail exist.
            Any qualifiers given will be used as a base for all SeqFeatures created.
        """
        # calculate locations
        get_sub = self.get_sub_location_from_protein_coordinates
        total_length = len(self.location) // 3
        if self._leader:
            leader_location = get_sub(0, len(self._leader))
        core_location = get_sub(len(self._leader), total_length - len(self._tail))
        if self._tail:
            tail_location = get_sub(total_length - len(self._tail), total_length)

        # build features
        features = []
        if self.leader:
            leader_qualifiers = {
                "prepeptide": ["leader"],
                "note": [
                    "peptide class: %s" % self.peptide_class,
                    "predicted sequence: %s" % self.leader,
                ],
                "locus_tag": [self.locus_tag],
                "aSTool": [self.tool],
                "tool": ["antismash"],
            }
            leader = SeqFeature(leader_location, type=self.type)
            leader.qualifiers.update(leader_qualifiers)
            leader.translation = self.leader
            features.append(leader)

        # use provided qualifiers, if exists
        if not qualifiers:
            qualifiers = {'note': []}
        if 'note' not in qualifiers:
            qualifiers['note'] = []
        core = SeqFeature(core_location, type=self.type, qualifiers=qualifiers)
        core.qualifiers.update({
            "prepeptide": ["core"],
            "locus_tag": [self.locus_tag],
            "peptide": [self.peptide_class],
            "predicted_class": [self.peptide_subclass],
            "score": ["{:.2f}".format(self.score)],
            "molecular_weight": ["{:.1f}".format(self.molecular_weight)],
            "monoisotopic_mass": ["{:.1f}".format(self.monoisotopic_mass)],
            "core_sequence": [self._core],
            "aSTool": [self.tool],
            "tool": ["antismash"],
        })
        if self._leader:
            core.qualifiers.update({
                "leader_sequence": [self._leader],
                "leader_location": [str(leader_location)],
            })
        if self._tail:
            core.qualifiers.update({
                "tail_sequence": [self._tail],
                "tail_location": [str(tail_location)],
            })
        if self.alternative_weights:
            weights = map(lambda x: "%0.1f" % x, self.alternative_weights)
            core.qualifiers['alternative_weights'] = list(weights)

        features.append(core)

        if self.tail:
            tail = SeqFeature(tail_location, type="CDS_motif")
            tail.translation = self.tail
            tail.qualifiers['locus_tag'] = [self.locus_tag]
            tail.qualifiers['note'] = ['tail peptide', self.peptide_class]
            tail.qualifiers['prepeptide'] = ['tail']
            tail.qualifiers['aSTool'] = [self.tool]
            tail.qualifiers["tool"] = ["antismash"]
            features.append(tail)

        return features

    @staticmethod
    def from_biopython(bio_feature: SeqFeature, feature: "Prepeptide" = None,  # type: ignore
                       leftovers: Dict[str, Any] = None) -> "Prepeptide":
        if leftovers is None:
            leftovers = Feature.make_qualifiers_copy(bio_feature)

        section = leftovers.pop("prepeptide", [""])[0]
        if not section:
            raise SecmetInvalidInputError("cannot reconstruct Prepeptide from biopython feature %s" % bio_feature)
        elif section != "core":
            raise SecmetInvalidInputError("Prepeptide can only be reconstructed from core feature")
        peptide_class = leftovers.pop("peptide")[0]
        core = leftovers.pop("core_sequence")[0]
        alt_weights = [float(weight) for weight in leftovers.pop("alternative_weights", [])]
        leader = leftovers.pop("leader_sequence", [""])[0]
        locations = [bio_feature.location]
        if leader:
            leader_location = location_from_string(leftovers.pop("leader_location")[0])
            locations.insert(0, leader_location)
        tail = leftovers.pop("tail_sequence", [""])[0]
        if tail:
            tail_location = location_from_string(leftovers.pop("tail_location")[0])
            locations.append(tail_location)

        location = build_location_from_others(locations)

        return Prepeptide(location, peptide_class, core,
                          leftovers.pop("locus_tag")[0],
                          leftovers.pop("aSTool")[0],
                          leftovers.pop("predicted_class")[0],
                          float(leftovers.pop("score")[0]),
                          float(leftovers.pop("monoisotopic_mass")[0]),
                          float(leftovers.pop("molecular_weight")[0]), alt_weights,
                          leader, tail)

    def to_json(self) -> Dict[str, Any]:
        """ Converts the qualifier to a dictionary for storing in JSON results.
        """
        data = dict(vars(self))
        for var in ["_tail", "_core", "_leader"]:
            data[var.replace("_", "")] = data[var]
            del data[var]
        data["location"] = str(self.location)
        data["score"] = self.score
        data["locus_tag"] = self.locus_tag
        data["tool"] = self.tool
        del data["detailed_information"]
        if self.detailed_information:
            data["detailed_info"] = self.detailed_information.to_biopython_qualifiers()
        return data

    @classmethod
    def from_json(cls, data: Dict[str, Any]) -> "Prepeptide":
        """ Rebuilds a Prepeptide instance from a JSON representation """
        details = data.pop("detailed_info", None)
        data["location"] = location_from_string(data["location"])
        peptide = cls(**data)
        peptide.detailed_information = rebuild_qualifier(details, peptide.peptide_class)
        return peptide
