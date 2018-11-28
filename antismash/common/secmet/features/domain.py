# License: GNU Affero General Public License v3 or later
# A copy of GNU AGPL v3 should have been included in this software package in LICENSE.txt.

""" A base class for all domain sub-features """

from collections import OrderedDict
from typing import Dict, List, Optional

from Bio.SeqFeature import SeqFeature

from antismash.common.secmet.qualifiers import ActiveSiteFinderQualifier

from .feature import Feature, Location
from .antismash_feature import AntismashFeature


class Domain(AntismashFeature):
    """ A base class for features which represent a domain type """
    __slots__ = ["domain", "_asf"]

    def __init__(self, location: Location, feature_type: str,
                 domain: Optional[str] = None, tool: str = None,
                 created_by_antismash: bool = True) -> None:
        super().__init__(location, feature_type, tool=tool, created_by_antismash=created_by_antismash)
        if domain is not None:
            if not isinstance(domain, str):
                raise TypeError("Domain must be given domain as a string, not %s" % type(domain))
            if not domain:
                raise ValueError("Domain cannot be an empty string")
        self.domain = domain
        self._asf = ActiveSiteFinderQualifier()

    @property
    def asf(self) -> ActiveSiteFinderQualifier:
        """ An ActiveSiteFinderQualifier storing active site descriptions """
        return self._asf

    def to_biopython(self, qualifiers: Dict[str, List[str]] = None) -> List[SeqFeature]:
        mine = OrderedDict()  # type: Dict[str, List[str]]
        if self.domain:
            mine["aSDomain"] = [self.domain]
        if self._asf:
            mine["ASF"] = self.asf.to_biopython()
        if qualifiers:
            mine.update(qualifiers)
        return super().to_biopython(mine)

    @staticmethod
    def from_biopython(bio_feature: SeqFeature, feature: "Domain" = None,  # type: ignore
                       leftovers: Dict[str, List[str]] = None) -> "Domain":
        if leftovers is None:
            leftovers = Feature.make_qualifiers_copy(bio_feature)
        if not feature:
            raise ValueError("Domain shouldn't be instantiated directly")
        else:
            assert isinstance(feature, Domain), type(feature)

        # grab optional qualifiers
        feature.domain = leftovers.pop("aSDomain", [""])[0] or None
        for asf_label in leftovers.pop("ASF", []):
            feature.asf.add(asf_label)

        # grab parent optional qualifiers
        updated = super(Domain, feature).from_biopython(bio_feature, feature=feature, leftovers=leftovers)
        assert updated is feature
        assert isinstance(updated, Domain)
        return updated
