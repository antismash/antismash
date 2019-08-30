# License: GNU Affero General Public License v3 or later
# A copy of GNU AGPL v3 should have been included in this software package in LICENSE.txt.

""" A feature to represent a Pfam domain """

from collections import OrderedDict
from typing import Any, Dict, List, Optional, Type, TypeVar

from Bio.SeqFeature import SeqFeature

from antismash.common.secmet.locations import FeatureLocation
from antismash.common.secmet.qualifiers import GOQualifier

from ..errors import SecmetInvalidInputError
from .feature import Feature, Location
from .domain import Domain

T = TypeVar("T", bound="PFAMDomain")


class PFAMDomain(Domain):
    """ A feature representing a PFAM domain within a CDS.
    """
    __slots__ = ["description", "probability",
                 "gene_ontologies", "identifier", "version"]
    FEATURE_TYPE = "PFAM_domain"

    def __init__(self, location: Location, description: str, protein_location: Location,
                 identifier: str, tool: str, locus_tag: str, domain: Optional[str] = None,
                 ) -> None:
        """ Arguments:
                location: the DNA location of the feature
                description: a string with a description
                protein_location: the location within the parent CDS translation
                identifier: the Pfam identifier (e.g. PF00067 or PF00067.14)
                tool: the name of the tool used to find/create the feature
                locus_tag: the name of the parent CDS feature
                domain: the name for the domain (e.g. p450 or 'Type III restriction enzyme')
        """
        super().__init__(location, self.FEATURE_TYPE, protein_location, locus_tag,
                         domain=domain, tool=tool)
        if not isinstance(description, str):
            raise TypeError("PFAMDomain description must be a string, not %s" % type(description))
        if not description:
            raise ValueError("PFAMDomain description cannot be empty")
        self.description = description
        self.probability = None  # type: Optional[float]
        self.version = None
        if not identifier:
            raise ValueError("Pfam identifier cannot be empty")
        if "." in identifier:
            identifier, version = identifier.split(".", maxsplit=1)
            self.version = int(version)
        if not (len(identifier) == 7 and identifier.startswith('PF') and identifier[2:].isdecimal()):
            raise ValueError("invalid Pfam identifier: %s" % identifier)
        self.identifier = str(identifier)
        self.gene_ontologies = None  # type: Optional[GOQualifier]

    @property
    def full_identifier(self) -> str:
        """ Returns the Pfam identifier with version, if available, in the form:
            PF00067.1
        """
        if not self.version:
            return self.identifier
        return "%s.%d" % (self.identifier, self.version)

    def to_biopython(self, qualifiers: Dict[str, List[str]] = None) -> List[SeqFeature]:
        mine = OrderedDict()  # type: Dict[str, List[str]]
        mine["description"] = [self.description]
        if self.probability is not None:
            mine["probability"] = [str(self.probability)]
        mine["db_xref"] = [self.full_identifier]
        if self.gene_ontologies:  # should only be the case if db_xrefs present, since those are needed for mapping
            mine["gene_ontologies"] = self.gene_ontologies.to_biopython()
            mine["db_xref"].extend(sorted(self.gene_ontologies.ids))
        if qualifiers:
            mine.update(qualifiers)
        return super().to_biopython(mine)

    @classmethod
    def from_biopython(cls: Type[T], bio_feature: SeqFeature, feature: T = None,
                       leftovers: Dict[str, List[str]] = None, record: Any = None) -> T:
        if leftovers is None:
            leftovers = Feature.make_qualifiers_copy(bio_feature)
        # grab mandatory qualifiers and create the class
        description = leftovers.pop("description")[0]
        p_start = int(leftovers.pop("protein_start")[0])
        p_end = int(leftovers.pop("protein_end")[0])
        xref = leftovers.get("db_xref", [])  # only remove the interesting part
        name = None
        for i, ref in enumerate(xref):
            if ref.startswith("PF"):
                name = ref
            xref.pop(i)
            break
        if name is None:
            raise SecmetInvalidInputError("PFAMDomain missing identifier")
        tool = leftovers.pop("aSTool")[0]
        locus_tag = leftovers.pop("locus_tag", ["(unknown)"])[0]

        feature = cls(bio_feature.location, description, FeatureLocation(p_start, p_end),
                      identifier=name, tool=tool, locus_tag=locus_tag)

        # grab optional qualifiers
        feature.gene_ontologies = GOQualifier.from_biopython(leftovers.pop("gene_ontologies", []))
        if "probability" in leftovers:
            feature.probability = float(leftovers["probability"][0])

        # grab parent optional qualifiers
        updated = super().from_biopython(bio_feature, feature=feature, leftovers=leftovers, record=record)
        assert isinstance(updated, PFAMDomain)
        return updated
