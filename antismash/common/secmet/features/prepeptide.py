# License: GNU Affero General Public License v3 or later
# A copy of GNU AGPL v3 should have been included in this software package in LICENSE.txt.

""" A class for representing RiPP prepeptides """

from typing import Dict, List

from Bio.SeqFeature import SeqFeature

from .cds_motif import CDSMotif
from .feature import FeatureLocation


class Prepeptide(CDSMotif):
    """ A class representing a prepeptide. Used for tracking a multi-feature
        construction with a leader, core and tail. To allow for multiple types
        of prepeptide (e.g. lanthi- or sacti-peptides), only the core must exist.
    """
    def __init__(self, location: FeatureLocation, peptide_class: str, core: str, locus_tag: str,
                 peptide_subclass: str = None, score: float = 0., monoisotopic_mass: float = 0.,
                 molecular_weight: float = 0., alternative_weights: List[float] = None,
                 leader: str = "", tail: str = "") -> None:
        """
            Arguments:
                peptide_class: the kind of prepeptide, e.g. 'lanthipeptide', 'thiopeptide'
                core: the sequence of the core
                locus_tag: the locus tag to use for the feature
                prepeptide_subclass: the subclass of the prepeptide, e.g. 'Type II'
                score: the prepeptide score
                monoisotopic_mass: the monoisotopic mass of the prepeptide
                molecular_weight: the molecular weight of the prepeptide
                alternative_weights: a list of alternative weights for the prepeptide
                leader: the sequence of the leader, if it exists
                tail: the sequence of the tail, if it exists
        """
        super().__init__(location)
        for arg in [peptide_class, core, leader, tail]:
            assert isinstance(arg, str), type(arg)
        self._leader = leader
        self._core = core
        self._tail = tail
        self.locus_tag = locus_tag
        self.type = "CDS_motif"
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
        return self.locus_tag

    def to_biopython(self, qualifiers: Dict[str, List] = None) -> List[SeqFeature]:
        """ Generates up to three SeqFeatures, depending if leader and tail exist.
            Any qualifiers given will be used as a base for all SeqFeatures created.
        """
        # calculate core location
        core_start = self.location.start
        core_end = self.location.end
        if self.leader:
            core_start += len(self.leader) * 3
        if self.tail:
            core_end -= len(self.tail) * 3
        core_location = FeatureLocation(core_start, core_end, self.location.strand)

        # add qualifiers
        if not qualifiers:
            qualifiers = {'note': []}
        if 'note' not in qualifiers:
            qualifiers['note'] = []

        # build features
        features = []
        if self.leader:
            start = self.location.start
            leader_location = FeatureLocation(start, core_location.start, self.location.strand)
            leader = SeqFeature(leader_location, type="CDS_motif", qualifiers={"note": []})
            leader.translation = self.leader
            leader.qualifiers['locus_tag'] = [self.locus_tag]
            leader.qualifiers['note'].extend(['leader peptide', self.peptide_class,
                                              'predicted leader seq: %s' % self.leader])
            features.append(leader)

        core = SeqFeature(core_location, type="CDS_motif", qualifiers=qualifiers)
        core.qualifiers['locus_tag'] = [self.locus_tag]
        core.qualifiers['note'].extend(['core peptide', self.peptide_class,
                                        'predicted class: %s' % self.peptide_subclass,
                                        "predicted core seq: %s" % self.core,
                                        "score: %0.2f" % self.score,
                                        "molecular weight: %0.1f" % self.molecular_weight,
                                        "monoisotopic mass: %0.1f" % self.monoisotopic_mass])
        if self.alternative_weights:
            weights = map(lambda x: "%0.1f" % x, self.alternative_weights)
            core.qualifiers['note'].append('alternative weights: %s' % "; ".join(weights))

        features.append(core)

        if self.tail:
            tail_location = FeatureLocation(core_location.end, self.location.end, self.location.strand)
            tail = SeqFeature(tail_location, type="CDS_motif")
            tail.translation = self.tail
            tail.qualifiers['locus_tag'] = [self.locus_tag]
            tail.qualifiers['note'] = ['tail peptide', self.peptide_class]
            features.append(tail)

        return features

    def to_json(self) -> Dict:
        """ Converts the qualifier to a dictionary for storing in JSON results.
        """
        data = dict(vars(self))
        for var in ["_tail", "_core", "_leader"]:
            data[var.replace("_", "")] = data[var]
            del data[var]
        data["location"] = str(self.location)
        data["score"] = self.score
        return data
