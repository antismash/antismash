# License: GNU Affero General Public License v3 or later
# A copy of GNU AGPL v3 should have been included in this software package in LICENSE.txt.

""" A feature to represent a gene """

from typing import Any, Dict, List, Optional

from Bio.SeqFeature import SeqFeature

from .feature import Feature, FeatureLocation


class Gene(Feature):
    """ A feature representing a Gene (more general than a CDS) """
    __slots__ = ["_pseudo", "locus_tag", "gene_name"]

    def __init__(self, location: FeatureLocation, locus_tag: Optional[str] = None,
                 gene_name: Optional[str] = None, pseudo_gene: bool = False,
                 created_by_antismash: bool = False, qualifiers: Optional[Dict[str, List[str]]] = None) -> None:
        super().__init__(location, feature_type="gene",
                         created_by_antismash=created_by_antismash)
        self.locus_tag = str(locus_tag) if locus_tag else None
        self.gene_name = str(gene_name) if gene_name else None
        if not self.locus_tag and not self.gene_name:
            raise ValueError("Gene instances must have a locus tag or name")
        self._pseudo = bool(pseudo_gene)
        if self._pseudo:
            assert not created_by_antismash, "pseudo genes can only come from input files"
        if qualifiers:
            assert isinstance(qualifiers, dict)
            self._qualifiers.update(qualifiers)

    def get_name(self) -> str:
        """ Returns the locus tag or gene name of the gene, in that order """
        if self.locus_tag is not None:
            return self.locus_tag
        if self.gene_name is not None:
            return self.gene_name
        raise ValueError("names removed after construction")

    def is_pseudo_gene(self) -> bool:
        """ Was the gene marked as a pseudo-gene """
        return self._pseudo

    def to_biopython(self, qualifiers: Dict[str, Any] = None) -> SeqFeature:
        """ Construct a matching SeqFeature for this Gene """
        if not qualifiers:
            qualifiers = {}
        if self.locus_tag:
            qualifiers["locus_tag"] = [self.locus_tag]
        if self.gene_name:
            qualifiers["gene"] = [self.gene_name]
        if self._pseudo:
            qualifiers["pseudo"] = []
        return super().to_biopython(qualifiers)

    @staticmethod
    def from_biopython(bio_feature: SeqFeature, feature: "Gene" = None,  # type: ignore
                       leftovers: Dict[str, List[str]] = None) -> "Gene":
        if leftovers is None:
            leftovers = Feature.make_qualifiers_copy(bio_feature)
        # grab mandatory qualifiers and create the class
        locus = leftovers.pop("locus_tag", [""])[0] or None
        name = leftovers.pop("gene", [""])[0] or None
        if not (locus or name):
            name = "gene%s_%s" % (bio_feature.location.start, bio_feature.location.end)
        pseudo = "pseudo" in leftovers
        if pseudo:
            leftovers.pop("pseudo")
        feature = Gene(bio_feature.location, locus_tag=locus, gene_name=name, pseudo_gene=pseudo)
        super(Gene, feature).from_biopython(bio_feature, feature=feature, leftovers=leftovers)
        return feature
