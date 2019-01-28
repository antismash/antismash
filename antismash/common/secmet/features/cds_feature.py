# License: GNU Affero General Public License v3 or later
# A copy of GNU AGPL v3 should have been included in this software package in LICENSE.txt.

""" A class for coding sequence (CDS) features """

from collections import OrderedDict
import logging
from typing import Any, Dict, List, Optional
from typing import Union  # comment hints  # pylint: disable=unused-import

from Bio.Data import IUPACData
from Bio.SeqFeature import SeqFeature

from antismash.common.secmet import features  # comment hints  # pylint:disable=unused-import
from antismash.common.secmet.qualifiers import (
    GeneFunction,
    GeneFunctionAnnotations,
    NRPSPKSQualifier,
    SecMetQualifier,
)

from ..errors import SecmetInvalidInputError
from ..locations import AfterPosition, BeforePosition, Location
from .feature import Feature


_VALID_TRANSLATION_CHARS = set(IUPACData.protein_letters).union(set("X"))


def _sanitise_id_value(name: Optional[str]) -> Optional[str]:
    """ Ensures a name doesn't contain characters that will break external programs"""
    if name is None:
        return None
    name = str(name)
    illegal_chars = set("!\"#$%&()*+,:; \r\n\t=>?@[]^`'{|}/ ")
    for char in set(name).intersection(illegal_chars):
        name = name.replace(char, "_")
    return name


def _is_valid_translation_length(translation: str, location: Location) -> bool:
    """ Returns True if the given translation is contained within an explicit
        location or the location has an ambiguous end
    """
    if len(translation) * 3 <= len(location):
        return True
    if location.strand == -1 and isinstance(location.start, BeforePosition):
        return True
    return location.strand != -1 and isinstance(location.end, AfterPosition)


class CDSFeature(Feature):
    """ A feature representing a single CDS/gene. """
    __slots__ = ["_translation", "protein_id", "locus_tag", "gene", "product",
                 "transl_table", "_sec_met", "_gene_functions",
                 "unique_id", "_nrps_pks", "motifs", "region"]

    def __init__(self, location: Location, translation: str, locus_tag: str = None,
                 protein_id: str = None, product: str = "", gene: str = None,
                 translation_table: int = 1) -> None:
        super().__init__(location, feature_type="CDS")
        if location.strand not in [1, -1]:
            raise ValueError("Strand must be 1 or -1 for a CDS, not %s" % location.strand)
        # mandatory
        self._gene_functions = GeneFunctionAnnotations()

        if not (protein_id or locus_tag or gene):
            raise ValueError("CDSFeature requires at least one of: gene, protein_id, locus_tag")
        # semi-optional
        self.protein_id = _sanitise_id_value(protein_id)
        self.locus_tag = _sanitise_id_value(locus_tag)
        self.gene = _sanitise_id_value(gene)
        self.translation = str(translation)

        # optional
        if not isinstance(product, str):
            raise TypeError("product must be a string, not %s" % type(product))
        self.product = product
        self.transl_table = int(translation_table)
        self._sec_met = SecMetQualifier()
        self._nrps_pks = NRPSPKSQualifier(self.location.strand)

        self.motifs = []  # type: List[features.CDSMotif]

        # runtime-only data
        self.region = None  # type: Optional[features.Region]
        self.unique_id = None  # type: Optional[str] # set only when added to a record

    @property
    def gene_functions(self) -> GeneFunctionAnnotations:
        """ All gene function annotations for the CDS """
        return self._gene_functions

    @property
    def gene_function(self) -> GeneFunction:
        """ The likely gene function of the CDS, as determined by all annotated
            gene functions.
        """
        return self._gene_functions.get_classification()

    @property
    def sec_met(self) -> Optional[SecMetQualifier]:
        """ The qualifier containing secondary metabolite information for the
            CDSFeature.
        """
        return self._sec_met

    @sec_met.setter
    def sec_met(self, sec_met: SecMetQualifier) -> None:
        if sec_met is not None and not isinstance(sec_met, SecMetQualifier):
            raise TypeError("CDSFeature.sec_met can only be set to an instance of SecMetQualifier")
        self._sec_met = sec_met

    @property
    def nrps_pks(self) -> NRPSPKSQualifier:
        """ The NRPSPKSQualifier of the feature """
        return self._nrps_pks

    @nrps_pks.setter
    def nrps_pks(self, qualifier: NRPSPKSQualifier) -> None:
        if qualifier is not None and not isinstance(qualifier, NRPSPKSQualifier):
            raise TypeError("CDSFeature.nrps_pks can only be set to an instance of NRPSPKSQualifier")
        self._nrps_pks = qualifier

    @property
    def translation(self) -> str:
        """ The translation of the CDS, as a string of amino acids """
        return self._translation

    @translation.setter
    def translation(self, translation: str) -> None:
        if not translation:
            raise ValueError("valid translation required")
        invalid = set(translation) - _VALID_TRANSLATION_CHARS
        if invalid:
            raise ValueError("invalid translation characters: %s" % invalid)
        if not _is_valid_translation_length(translation, self.location):
            raise ValueError("translation longer than location allows: %s > %s" % (
                                len(translation) * 3, len(self.location)))
        self._translation = translation  # pylint: disable=attribute-defined-outside-init

    def get_accession(self) -> str:
        "Get the gene ID from protein id, gene name or locus_tag, in that order"
        for val in [self.protein_id, self.gene, self.locus_tag]:
            if val:
                return val
        raise ValueError("%s altered to contain no identifiers" % self)

    def get_name(self) -> str:
        "Get the gene ID from locus_tag, gene name or protein id, in that order"
        for val in [self.locus_tag, self.gene, self.protein_id]:
            if val:
                return val
        raise ValueError("%s altered to contain no identifiers" % self)

    @staticmethod
    def from_biopython(bio_feature: SeqFeature, feature: "CDSFeature" = None,  # type: ignore
                       leftovers: Optional[Dict] = None, record: Any = None) -> "CDSFeature":
        if leftovers is None:
            leftovers = Feature.make_qualifiers_copy(bio_feature)
        # grab mandatory qualifiers
        transl_table = 1
        if record:
            transl_table = record.transl_table
        if "transl_table" in leftovers:
            transl_table = int(leftovers.pop("transl_table")[0])
        translation = leftovers.pop("translation", [""])[0]

        # semi-optional qualifiers
        protein_id = leftovers.pop("protein_id", [None])[0]
        locus_tag = leftovers.pop("locus_tag", [None])[0]
        gene = leftovers.pop("gene", [None])[0]
        if not (gene or protein_id or locus_tag):
            if "pseudo" in leftovers or "pseudogene" in leftovers:
                gene = "pseudo%s_%s"
            else:
                gene = "cds%s_%s"
            gene = gene % (bio_feature.location.start, bio_feature.location.end)

        # ensure translation is valid if it exists
        if translation:
            invalid = set(translation) - _VALID_TRANSLATION_CHARS
            if invalid:
                logging.warning("Regenerating translation for CDS %s (at %s) containing invalid characters: %s",
                                locus_tag or protein_id or gene, bio_feature.location, invalid)
                translation = ""
        # ensure that the translation fits
        if not _is_valid_translation_length(translation, bio_feature.location):
            raise SecmetInvalidInputError("translation longer than location allows: %s > %s" % (
                                len(translation) * 3, len(bio_feature.location)))
        # finally, generate the translation if it doesn't exist
        if not translation:
            if not record:
                raise SecmetInvalidInputError("no translation in CDS and no record to generate it with")
            if bio_feature.location.end > len(record.seq):
                raise SecmetInvalidInputError("feature missing translation and sequence too short: %s" % (
                                              (gene or protein_id or locus_tag)))
            translation = record.get_aa_translation_from_location(bio_feature.location, transl_table)

        assert _is_valid_translation_length(translation, bio_feature.location)

        feature = CDSFeature(bio_feature.location, translation, gene=gene,
                             locus_tag=locus_tag, protein_id=protein_id,
                             translation_table=transl_table)

        # grab optional qualifiers
        feature.product = leftovers.pop("product", [""])[0]
        sec_met = leftovers.pop("sec_met_domain", None)
        if sec_met:
            feature.sec_met = SecMetQualifier.from_biopython(sec_met)
        gene_functions = leftovers.pop("gene_functions", [])
        if gene_functions:
            feature.gene_functions.add_from_qualifier(gene_functions)

        # grab parent optional qualifiers
        super(CDSFeature, feature).from_biopython(bio_feature, feature=feature, leftovers=leftovers)

        return feature

    def to_biopython(self, qualifiers: Dict[str, List[str]] = None) -> SeqFeature:
        mine = OrderedDict()  # type: Dict[str, List[str]]
        # mandatory
        mine["translation"] = [self.translation]
        # optional
        for attr in ["gene", "transl_table", "locus_tag",
                     "protein_id", "product"]:
            val = getattr(self, attr)
            if val:
                mine[attr] = [str(val)]
        if self._gene_functions:
            mine["gene_functions"] = list(map(str, self._gene_functions))
            mine["gene_kind"] = [str(self.gene_function)]
        if self.sec_met:
            mine["sec_met_domain"] = list(map(str, self.sec_met))
        # respect qualifiers given to us
        if qualifiers:
            mine.update(qualifiers)
        return super().to_biopython(mine)

    def __repr__(self) -> str:
        return str(self)

    def __str__(self) -> str:
        return "CDS(%s, %s)" % (self.get_name(), self.location)
