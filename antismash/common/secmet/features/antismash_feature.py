# License: GNU Affero General Public License v3 or later
# A copy of GNU AGPL v3 should have been included in this software package in LICENSE.txt.

""" A base class for all antiSMASH-specific features """

from collections import OrderedDict
from typing import Dict, List

from Bio.SeqFeature import SeqFeature

from .feature import Feature, FeatureLocation


class AntismashFeature(Feature):
    """ A base class for all sub-CDS Antismash features """
    __slots__ = ["domain_id", "database", "detection", "_evalue", "label",
                 "locus_tag", "_score", "_translation"]

    def __init__(self, location: FeatureLocation, feature_type: str) -> None:
        super().__init__(location, feature_type, created_by_antismash=True)
        self.domain_id = None  # type: str
        self.database = None  # type: str
        self.detection = None  # type: str
        self._evalue = None  # type: float
        self.label = None  # type: str
        self.locus_tag = None  # type: str
        self._score = None  # type: float

        self._translation = ""

    @property
    def translation(self) -> str:
        """ The amino acid translation of the feature. """
        if not self._translation:
            raise ValueError("Domain has no translation: %s" % self.domain_id)
        return self._translation

    @translation.setter
    def translation(self, translation: str) -> None:
        assert isinstance(translation, str)
        if not translation:
            raise ValueError("Domain translation cannot be empty")
        self._translation = translation

    @property
    def score(self) -> float:
        """ The bitscore reported by a tool when locating the feature """
        return self._score

    @score.setter
    def score(self, score: float) -> None:
        self._score = float(score)

    @property
    def evalue(self) -> float:
        """ The e-value reported by a tool when locating the feature """
        return self._evalue

    @evalue.setter
    def evalue(self, evalue: float) -> None:
        self._evalue = float(evalue)

    def get_name(self) -> str:
        """ Returns the domain's identifier """
        return self.domain_id

    def to_biopython(self, qualifiers: Dict[str, List[str]] = None) -> List[SeqFeature]:
        mine = OrderedDict()  # type: Dict[str, List[str]]
        if self.label:
            mine["label"] = [self.label]
        if self.score is not None:
            mine["score"] = [str(self.score)]
        if self.evalue is not None:
            mine["evalue"] = [str("{:.2E}".format(self.evalue))]
        if self.locus_tag:
            mine["locus_tag"] = [self.locus_tag]
        if self._translation:
            mine["translation"] = [self._translation]
        if self.database:
            mine["database"] = [self.database]
        if self.detection:
            mine["detection"] = [self.detection]
        if self.domain_id:
            mine["domain_id"] = [self.domain_id]
        if qualifiers:
            mine.update(qualifiers)
        return super().to_biopython(mine)

    @staticmethod
    def from_biopython(bio_feature: SeqFeature, feature: "AntismashFeature" = None,  # type: ignore
                       leftovers: Dict[str, List[str]] = None) -> "AntismashFeature":
        if leftovers is None:
            leftovers = Feature.make_qualifiers_copy(bio_feature)
        if not feature:
            raise ValueError("AntismashFeature shouldn't be instantiated directly")
        else:
            assert isinstance(feature, AntismashFeature)
        # grab optional qualifiers
        feature.domain_id = leftovers.pop("domain_id", [None])[0]
        feature.database = leftovers.pop("database", [None])[0]
        feature.detection = leftovers.pop("detection", [None])[0]
        feature.label = leftovers.pop("label", [None])[0]
        feature.locus_tag = leftovers.pop("locus_tag", [None])[0]
        translation = leftovers.pop("translation", [None])[0]
        if translation is not None:
            feature.translation = translation
        if "evalue" in leftovers:
            feature.evalue = float(leftovers.pop("evalue")[0])
        if "score" in leftovers:
            feature.score = float(leftovers.pop("score")[0])

        # grab parent optional qualifiers
        updated = super(AntismashFeature, feature).from_biopython(bio_feature, feature=feature, leftovers=leftovers)
        assert isinstance(updated, AntismashFeature)
        return updated
