# License: GNU Affero General Public License v3 or later
# A copy of GNU AGPL v3 should have been included in this software package in LICENSE.txt.

""" A feature to represent a Pfam domain """

from collections import OrderedDict
from typing import Dict, List, Optional

from Bio.SeqFeature import SeqFeature

from antismash.common.secmet.qualifiers import GOQualifier

from .feature import Feature, FeatureLocation
from .domain import Domain


class PFAMDomain(Domain):
    """ A feature representing a PFAM domain within a CDS.
    """
    __slots__ = ["description", "db_xref", "probability", "protein_start", "protein_end",
                 "gene_ontologies"]

    def __init__(self, location: FeatureLocation, description: str, protein_start: int,
                 protein_end: int, domain: Optional[str] = None) -> None:
        """ Arguments:
                location: the DNA location of the feature
                description: a string with a description
                protein_start: the start point within the parent CDS translation
                protein_end: the end point within the parent CDS translation
                domain: the name for the domain (e.g. p450 vs the dbxref PF00067)
        """
        super().__init__(location, feature_type="PFAM_domain", domain=domain)
        if not isinstance(description, str):
            raise TypeError("PFAMDomain description must be a string, not %s" % type(description))
        if not description:
            raise ValueError("PFAMDomain description cannot be empty")
        self.description = description
        self.probability = None
        self.db_xref = []  # type: List[str]
        self.protein_start = int(protein_start)
        self.protein_end = int(protein_end)
        if self.protein_start >= self.protein_end:
            raise ValueError("A PFAMDomain protein location cannot end before it starts")
        self.gene_ontologies = None  # type: Optional[GOQualifier]

    def to_biopython(self, qualifiers: Dict[str, List[str]] = None) -> List[SeqFeature]:
        mine = OrderedDict()  # type: Dict[str, List[str]]
        mine["description"] = [self.description]
        mine["protein_start"] = [str(self.protein_start)]
        mine["protein_end"] = [str(self.protein_end)]
        if self.probability is not None:
            mine["probability"] = [self.probability]
        if self.db_xref:
            mine["db_xref"] = self.db_xref
        if self.gene_ontologies:  # should only be the case if db_xrefs present, since those are needed for mapping
            mine["gene_ontologies"] = self.gene_ontologies.to_biopython()
            mine["db_xref"].extend(sorted(self.gene_ontologies.ids))
        if qualifiers:
            mine.update(qualifiers)
        return super().to_biopython(mine)

    @staticmethod
    def from_biopython(bio_feature: SeqFeature, feature: "PFAMDomain" = None,  # type: ignore
                       leftovers: Dict[str, List[str]] = None) -> "PFAMDomain":
        if leftovers is None:
            leftovers = Feature.make_qualifiers_copy(bio_feature)
        # grab mandatory qualifiers and create the class
        description = leftovers.pop("description")[0]
        p_start = int(leftovers.pop("protein_start")[0])
        p_end = int(leftovers.pop("protein_end")[0])
        feature = PFAMDomain(bio_feature.location, description, p_start, p_end)

        # grab optional qualifiers
        feature.db_xref = leftovers.pop("db_xref", [])
        feature.gene_ontologies = GOQualifier.from_biopython(leftovers.pop("gene_ontologies", []))

        # grab parent optional qualifiers
        updated = super(PFAMDomain, feature).from_biopython(bio_feature, feature=feature, leftovers=leftovers)
        assert isinstance(updated, PFAMDomain)
        return updated
