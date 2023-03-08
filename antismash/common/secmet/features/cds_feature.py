# License: GNU Affero General Public License v3 or later
# A copy of GNU AGPL v3 should have been included in this software package in LICENSE.txt.

""" A class for coding sequence (CDS) features """

from collections import OrderedDict
import logging
from typing import Any, Dict, List, Optional, Tuple, Type, TypeVar

from Bio.Data import CodonTable, IUPACData
from Bio.SeqFeature import SeqFeature

from antismash.common.secmet import features
from antismash.common.secmet.qualifiers import (
    GeneFunction,
    GeneFunctionAnnotations,
    NRPSPKSQualifier,
    SecMetQualifier,
)

from ..errors import SecmetInvalidInputError
from ..locations import (
    AfterPosition,
    BeforePosition,
    frameshift_location_by_qualifier,
    Location,
)
from .feature import Feature, pop_locus_qualifier
from .module import Module

_VALID_TRANSLATION_CHARS = set(IUPACData.extended_protein_letters)
T = TypeVar("T", bound="CDSFeature")

MAX_TRANSLATION_LENGTH = 100_000


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


def _translation_fits_in_record(translation_length: int, location: Location,
                                record_length: int) -> bool:
    """ Checks that a translation fits within a record based on the feature location
        Allows for a single ambiguous amino to overlap the edges of a record

        Arguments:
            translation_length: length of translation in nucleotides
            location: the Location of the feature
            record_length: the length of the record in nucleotides

        Returns:
            True if location and translation length fit into the given record length
    """
    extra_length = translation_length - len(location)
    if extra_length <= 0:
        return True  # assumes the feature is contained by the record, of course

    if location.strand == -1:
        return location.start - extra_length >= -2

    return extra_length + location.end <= record_length + 2


def _ensure_valid_translation(translation: str, location: Location, transl_table: int,
                              record: Optional[Any]) -> str:
    """ Ensures that a given translation is valid for the matching location
        and record (if given). If a record is given and a translation contains
        invalid characters, an attempt will be made to generate a valid translation.

        Arguments:
            translation: the existing translation, if any
            location: the location of the feature owning the translation
            record: the Record the feature belongs to, or None to skip those checks and regeneration

        Returns:
            a valid translation
    """
    # ensure translation is valid if it exists
    if translation:
        invalid = set(translation) - _VALID_TRANSLATION_CHARS
        if invalid:
            logging.warning("Regenerating translation for CDS at %s containing invalid characters: %s",
                            location, invalid)
            translation = ""
    # ensure that the translation fits
    if not _is_valid_translation_length(translation, location):
        raise ValueError(f"translation longer than location allows: {len(translation) * 3} > {len(location)}")
    # if an arbitrary section of a record is used, the record can be too short for a given translation
    if record and not _translation_fits_in_record(len(translation)*3, location, len(record.seq)):
        raise ValueError("feature translation extends out of record")
    # finally, generate the translation if it doesn't exist
    if not translation:
        if not record:
            raise ValueError("no translation in CDS and no record to generate it with")
        if location.end > len(record.seq):
            raise ValueError("feature missing translation and sequence too short")
        if len(location) < 3:
            raise ValueError("CDS too short to generate translation")
        try:
            translation = record.get_aa_translation_from_location(location, transl_table)
        except CodonTable.TranslationError as err:
            raise ValueError(f"invalid codon: {err}")

    if len(translation) >= MAX_TRANSLATION_LENGTH:
        raise ValueError(f"translation too long for dependencies: {len(translation)}")

    assert _is_valid_translation_length(translation, location)
    return translation


def _verify_location(location: Location) -> None:
    """ Raises a relevant exception if the location is invalid for a CDS """
    if location.strand not in [1, -1]:
        if len(location.parts) > 1 and location.parts[0].strand in [1, -1]:
            raise ValueError("compound locations with mixed strands are not supported")
        raise ValueError(f"invalid strand: {location.strand}")


class CDSFeature(Feature):
    """ A feature representing a single CDS/gene. """
    __slots__ = ["_translation", "protein_id", "locus_tag", "gene", "product",
                 "transl_table", "_sec_met", "_gene_functions", "_modules",
                 "unique_id", "_nrps_pks", "motifs", "region"]
    FEATURE_TYPE = "CDS"

    def __init__(self, location: Location, translation: str, locus_tag: str = None,
                 protein_id: str = None, product: str = "", gene: str = None,
                 translation_table: int = 1) -> None:
        super().__init__(location, feature_type=self.FEATURE_TYPE)
        _verify_location(location)
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
            raise TypeError(f"product must be a string, not {type(product)}")
        self.product = product
        self.transl_table = int(translation_table)
        self._sec_met = SecMetQualifier()
        self._nrps_pks = NRPSPKSQualifier(self.location.strand)

        self._modules: List[Module] = []
        self.motifs: List[features.CDSMotif] = []

        # runtime-only data
        self.region: Optional[features.Region] = None
        self.unique_id: Optional[str] = None  # set only when added to a record

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
    def sec_met(self) -> SecMetQualifier:
        """ The qualifier containing secondary metabolite information for the
            CDSFeature.
        """
        return self._sec_met

    @sec_met.setter
    def sec_met(self, sec_met: SecMetQualifier) -> None:
        if not isinstance(sec_met, SecMetQualifier):
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
        if len(translation) >= MAX_TRANSLATION_LENGTH:
            detail = f"'{self.get_name()}': {len(translation)} aminos"
            raise ValueError(f"translation too long for dependencies in {detail}")
        invalid = set(translation) - _VALID_TRANSLATION_CHARS
        if invalid:
            raise ValueError(f"invalid translation characters: {invalid}")
        if not _is_valid_translation_length(translation, self.location):
            raise ValueError(f"translation longer than location allows: {len(translation) * 3} > {len(self.location)}")
        self._translation = translation  # pylint: disable=attribute-defined-outside-init

    @property
    def modules(self) -> Tuple[Module, ...]:
        """ Returns any modules present in the CDS. Modules may be included in
            more than one CDS feature.
        """
        return tuple(self._modules)

    def add_module(self, module: Module) -> None:
        """ Add a Module feature to the CDS, the module must overlap with the CDS """
        assert module.overlaps_with(self)
        self._modules.append(module)

    def get_accession(self) -> str:
        "Get the gene ID from protein id, gene name or locus_tag, in that order"
        for val in [self.protein_id, self.gene, self.locus_tag]:
            if val:
                return val
        raise ValueError(f"{self} altered to contain no identifiers")

    def get_name(self) -> str:
        "Get the gene ID from locus_tag, gene name or protein id, in that order"
        for val in [self.locus_tag, self.gene, self.protein_id]:
            if val:
                return val
        raise ValueError(f"{self} altered to contain no identifiers")

    @classmethod
    def from_biopython(cls: Type[T], bio_feature: SeqFeature, feature: T = None,
                       leftovers: Optional[Dict] = None, record: Any = None) -> T:
        if leftovers is None:
            leftovers = Feature.make_qualifiers_copy(bio_feature)
        # grab mandatory qualifiers
        transl_table = 1
        if record:
            transl_table = record.transl_table

        # semi-optional qualifiers
        protein_id = leftovers.pop("protein_id", [None])[0]
        locus_tag = pop_locus_qualifier(leftovers, allow_missing=True, default=None)
        gene = leftovers.pop("gene", [None])[0]
        if not (gene or protein_id or locus_tag):
            if "pseudo" in leftovers or "pseudogene" in leftovers:
                gene = "pseudo%s_%s"
            else:
                gene = "cds%s_%s"
            gene = gene % (bio_feature.location.start, bio_feature.location.end)
        name = locus_tag or protein_id or gene

        if "transl_table" in leftovers:
            raw_table = leftovers.pop("transl_table")[0]
            try:
                transl_table = int(raw_table)
            except ValueError:
                raise SecmetInvalidInputError(f"invalid translation table {raw_table!r} for CDS {name!r}")

        try:
            _verify_location(bio_feature.location)
        except Exception as err:
            message = f"invalid location for {name}: {err}"
            raise SecmetInvalidInputError(message) from err

        try:
            # before extracting a new translation, ensure that the location used
            # is correctly adjusted to account for a "codon_start" qualifier
            location = bio_feature.location
            if "codon_start" in leftovers:
                location = frameshift_location_by_qualifier(location, leftovers["codon_start"][0])
            translation = _ensure_valid_translation(leftovers.pop("translation", [""])[0],
                                                    location, transl_table, record)
        except ValueError as err:
            raise SecmetInvalidInputError(f"{err}: {name}") from err

        feature = cls(bio_feature.location, translation, gene=gene,
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
        feature.nrps_pks.add_from_qualifier(leftovers.pop("NRPS_PKS", []))

        # grab parent optional qualifiers
        super().from_biopython(bio_feature, feature=feature, leftovers=leftovers, record=record)

        return feature

    def to_biopython(self, qualifiers: Dict[str, List[str]] = None) -> SeqFeature:
        mine: Dict[str, List[str]] = OrderedDict()
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
        if self.nrps_pks:
            mine["NRPS_PKS"] = list(map(str, self.nrps_pks))
        # respect qualifiers given to us
        if qualifiers:
            mine.update(qualifiers)
        return super().to_biopython(mine)

    def __repr__(self) -> str:
        return str(self)

    def __str__(self) -> str:
        return f"CDS({self.get_name()}, {self.location})"

    def strip_antismash_annotations(self) -> None:
        """ Remove all antiSMASH-specific annotations from the feature """
        self.sec_met = SecMetQualifier()
        self.gene_functions.clear()
        self.nrps_pks = NRPSPKSQualifier(self.location.strand)
        self._modules.clear()
