# License: GNU Affero General Public License v3 or later
# A copy of GNU AGPL v3 should have been included in this software package in LICENSE.txt.

""" A base class for all domain sub-features """

from collections import OrderedDict
from typing import Any, Dict, List, Optional, Type, TypeVar

from Bio.SeqFeature import SeqFeature

from antismash.common.secmet.qualifiers import ActiveSiteFinderQualifier

from .feature import Feature, FeatureLocation, Location
from .antismash_feature import AntismashFeature

T = TypeVar("T", bound="Domain")


def generate_protein_location_from_qualifiers(qualifiers: Dict[str, List[str]], record: Any) -> FeatureLocation:
    """ Generates a protein location from the expected qualifiers of a Domain subclass.
        If older qualifiers are provided and the containing record is provided,
        including the CDS containing the domain, the location will be generated.

        Arguments:
            qualifiers: the raw qualifiers of the feature, expecting:
                     "protein_start" and "protein_end" for the simple case
                     "locus_tag" and "translation" for the regeneration case
            record: optionally, the containing Record

        Returns:
            the location within the parent CDS translation as a FeatureLocation
    """
    raw_start = qualifiers.get("protein_start", [""])[0]
    raw_end = qualifiers.get("protein_end", [""])[0]

    if raw_start:
        start = int(raw_start)
    else:
        if not record:
            raise ValueError("missing record for protein location regeneration")
        parent = record.get_cds_by_name(qualifiers["locus_tag"][0])
        start = parent.translation.find(qualifiers["translation"][0])
        raw_end = ""  # don't trust an end with no start

    if raw_end:
        end = int(raw_end)
    else:
        end = start + len(qualifiers["translation"])

    return FeatureLocation(start, end)


class Domain(AntismashFeature):
    """ A base class for features which represent a domain type """
    __slots__ = ["domain", "_asf", "protein_location"]

    def __init__(self, location: Location, feature_type: str, protein_location: FeatureLocation,
                 locus_tag: str, domain: Optional[str] = None, tool: str = None,
                 created_by_antismash: bool = True) -> None:
        super().__init__(location, feature_type, tool=tool, created_by_antismash=created_by_antismash)
        if domain is not None:
            if not isinstance(domain, str):
                raise TypeError("Domain must be given domain as a string, not %s" % type(domain))
            if not domain:
                raise ValueError("Domain cannot be an empty string")
        self.domain = domain
        self.protein_location = protein_location
        if not isinstance(protein_location, FeatureLocation):
            raise TypeError("protein location must be a FeatureLocation, not %s" % type(protein_location))
        if not locus_tag:
            raise ValueError("locus tag cannot be an empty string")
        self.locus_tag: str = str(locus_tag)
        self._asf = ActiveSiteFinderQualifier()

    @property
    def asf(self) -> ActiveSiteFinderQualifier:
        """ An ActiveSiteFinderQualifier storing active site descriptions """
        return self._asf

    def to_biopython(self, qualifiers: Dict[str, List[str]] = None) -> List[SeqFeature]:
        mine: Dict[str, List[str]] = OrderedDict()
        mine["protein_start"] = [str(int(self.protein_location.start))]
        mine["protein_end"] = [str(int(self.protein_location.end))]
        if self.domain:
            mine["aSDomain"] = [self.domain]
        if self._asf:
            mine["ASF"] = self.asf.to_biopython()
        if qualifiers:
            mine.update(qualifiers)
        return super().to_biopython(mine)

    @classmethod
    def from_biopython(cls: Type[T], bio_feature: SeqFeature, feature: T = None,
                       leftovers: Dict[str, List[str]] = None, record: Any = None) -> T:
        if leftovers is None:
            leftovers = Feature.make_qualifiers_copy(bio_feature)
        if not feature:
            raise ValueError("Domain shouldn't be instantiated directly")
        assert isinstance(feature, Domain), type(feature)

        # clean up qualifiers that must have been used already
        leftovers.pop("protein_start", None)
        leftovers.pop("protein_end", None)

        # grab optional qualifiers
        feature.domain = leftovers.pop("aSDomain", [""])[0] or None
        for asf_label in leftovers.pop("ASF", []):
            feature.asf.add(asf_label)

        # grab parent optional qualifiers
        updated = super().from_biopython(bio_feature, feature=feature, leftovers=leftovers, record=record)
        assert updated is feature
        assert isinstance(updated, Domain)
        return updated
