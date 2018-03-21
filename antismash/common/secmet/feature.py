# License: GNU Affero General Public License v3 or later
# A copy of GNU AGPL v3 should have been included in this software package in LICENSE.txt.

""" A collection of classes for representing a variety of feature types """

from collections import OrderedDict
import logging
import os
import warnings
from typing import Any, Dict, Iterable, List, Optional, Tuple, Union

from helperlibs.bio import seqio

from Bio.Seq import Seq
from Bio.SeqFeature import SeqFeature, FeatureLocation, CompoundLocation
from Bio.SeqRecord import SeqRecord

from .qualifiers import NRPSPKSQualifier, SecMetQualifier, GeneFunction, GeneFunctionAnnotations



def _convert_protein_position_to_dna(start: int, end: int, location: FeatureLocation) -> Tuple[int, int]:
    """ Convert a protein position to a nucleotide sequence position for use in generating
        new FeatureLocations from existing FeatureLocations and/or CompoundLocations.

        Arguments:
            position: the position in question, must be contained by the location
            location: the location of the related feature, for handling introns/split locations

        Returns:
            an int representing the calculated DNA location
    """
    if not 0 <= start < end <= len(location) // 3:
        raise ValueError("Protein positions %d and %d must be contained by %s" % (start, end, location))
    if location.strand == -1:
        dna_start = location.start + len(location) - end * 3
        dna_end = location.start + len(location) - start * 3
    else:
        dna_start = location.start + start * 3
        dna_end = location.start + end * 3

    # only CompoundLocations are complicated
    if not isinstance(location, CompoundLocation):
        assert location.start <= dna_start < dna_end <= location.end, "Converted coordinates %d..%d out of bounds for location %s" % (dna_start, dna_end, location)
        return dna_start, dna_end

    parts = sorted(location.parts, key=lambda x: x.start)
    gap = 0
    last_end = parts[0].start
    start_found = False
    end_found = False
    for part in parts:
        if start_found and end_found:
            break
        gap += part.start - last_end
        if not start_found and dna_start + gap in part:
            start_found = True
            dna_start = dna_start + gap
        if not end_found and dna_end + gap - 1 in part:
            end_found = True
            dna_end = dna_end + gap

        last_end = part.end

    assert start_found
    assert end_found

    assert location.start <= dna_start < dna_end <= location.end, "Converted coordinates %d..%d out of bounds for location %s" % (dna_start, dna_end, location)
    return dna_start, dna_end


class Feature:
    """ The base class of any feature. Contains only a location, the label of the
        subclass, the 'notes' qualifier, and other qualifiers not tracked by any
        subclass.
    """
    __slots__ = ["location", "notes", "type", "_qualifiers", "created_by_antismash"]

    def __init__(self, location: FeatureLocation, feature_type: str,
                 created_by_antismash: bool = False) -> None:
        assert isinstance(location, (FeatureLocation, CompoundLocation)), type(location)
        assert location.start < location.end, "Feature location invalid"
        self.location = location
        self.notes = []  # type: List[str]
        assert feature_type
        self.type = str(feature_type)
        self._qualifiers = OrderedDict()  # type: Dict[str, List[str]]
        self.created_by_antismash = bool(created_by_antismash)

    @property
    def strand(self) -> int:
        """ A simple wrapper to access the location strand directly.
        """
        return self.location.strand

    def extract(self, sequence: Seq) -> Seq:
        """ Extracts the section of the given sequence that this feature covers.

            Return type is always a Seq, unlike location.extract.
        """
        assert isinstance(sequence, Seq)
        return self.location.extract(sequence)

    def get_sub_location_from_protein_coordinates(self, start: int, end: int) -> FeatureLocation:
        """ Generates a FeatureLocation for a protein sequence based on the start
            and end positions within the features protein sequence.

            The start position is inclusive and the end position is exclusive.

            Arguments:
                start: a position between 0 and len(feature.location) // 3 - 1, inclusive
                end: a position between 1 and len(feature.location) // 3, inclusive
        """
        if not 0 <= start <= len(self.location) // 3 - 1:
            raise ValueError("Protein start coordinate must be contained by the feature")

        if not 1 <= end <= len(self.location) // 3:
            raise ValueError("Protein end coordinate must be contained by the feature")

        if start >= end:
            raise ValueError("Protein start coordinate must be less than the end coordinate")

        dna_start, dna_end = _convert_protein_position_to_dna(start, end, self.location)

        if not 0 <= dna_start - self.location.start < self.location.end - 2:
            raise ValueError("Protein coordinate start %d (nucl %d) is outside feature %s" % (start, dna_start, self))
        if not 2 < dna_end - self.location.start - 1 <= self.location.end:
            raise ValueError("Protein coordinate end %d (nucl %d) is outside feature %s" % (end, dna_end, self))

        if not isinstance(self.location, CompoundLocation):
            return FeatureLocation(dna_start, dna_end, self.location.strand)

        new_locations = []
        for location in sorted(self.location.parts, key=lambda x: x.start):
            if dna_start in location:
                new = FeatureLocation(dna_start, location.end, self.location.strand)
                # the end could also be in this part
                if dna_end - 1 in location:
                    # can't be a compound location with only one, so return a simple one
                    return FeatureLocation(new.start, dna_end, new.strand)
                new_locations.append(new)
            elif dna_end - 1 in location:  # 'in' uses start <= value < end
                new = FeatureLocation(location.start, dna_end, self.location.strand)
                new_locations.append(new)
                break
            elif new_locations:  # found a start, but haven't yet found an end
                new_locations.append(location)

        if not new_locations:
            raise ValueError("Could not create compound location from %s and internal protein coordinates %d..%d (dna %d..%d)" % (
                                str(self.location), start, end, dna_start, dna_end))
        if self.location.strand == -1:
            new_locations.reverse()

        if len(new_locations) == 1:
            return new_locations[0]

        return CompoundLocation(new_locations)

    def get_qualifier(self, key: str) -> Optional[Tuple]:
        """ Fetches a qualifier by key and returns a tuple of items stored under
            that key or None if the key was not present.
        """
        qualifier = self._qualifiers.get(key)
        if qualifier:
            return tuple(qualifier)
        return None

    def overlaps_with(self, other: Union["Feature", FeatureLocation]) -> bool:
        """ Returns True if the given feature overlaps with this feature.
            This operation is commutative, a.overlaps_with(b) is equivalent to
            b.overlaps_with(a).
        """
        if isinstance(other, Feature):
            location = other.location
        elif isinstance(other, FeatureLocation):
            location = other
        else:
            raise TypeError("Container must be a Feature or a FeatureLocation, not %s" % type(other))
        return (self.location.start in location
                or self.location.end - 1 in location
                or location.start in self.location
                or location.end - 1 in self.location)

    def is_contained_by(self, other: Union["Feature", FeatureLocation]) -> bool:
        """ Returns True if the given feature is wholly contained by this
            feature.
        """
        end = self.location.end - 1  # to account for the non-inclusive end
        if isinstance(other, Feature):
            return self.location.start in other.location and end in other.location
        if isinstance(other, FeatureLocation):
            return self.location.start in other and end in other
        raise TypeError("Container must be a Feature or a FeatureLocation, not %s" % type(other))

    def to_biopython(self, qualifiers: Dict[str, Any] = None) -> List[SeqFeature]:
        """ Converts this feature into one or more SeqFeature instances.

            Subclasses must manage their own attributes and potential extra
            features.
        """
        feature = SeqFeature(self.location, type=self.type)
        quals = self._qualifiers.copy()
        notes = self._qualifiers.get("note", []) + self.notes
        if qualifiers:
            notes += qualifiers.pop("note", [])
            quals.update(qualifiers)
        if notes:
            # sorting helps with consistency and comparison testing
            quals["note"] = sorted(notes)
        if self.created_by_antismash:
            quals["tool"] = ["antismash"]
        # sorted here to match the behaviour of biopython
        for key, val in sorted(quals.items()):
            feature.qualifiers[key] = val
        assert isinstance(feature.qualifiers, dict)
        return [feature]

    def __lt__(self, other: "Feature") -> bool:
        """ Allows sorting Features by location without key complication """
        assert isinstance(other, Feature)
        if self.location.start < other.location.start:
            return True
        elif self.location.start == other.location.start:
            return self.location.end < other.location.end
        return False

    def __str__(self) -> str:
        return repr(self)

    def __repr__(self) -> str:
        return "%s(%s)" % (self.type, self.location)

    @staticmethod
    def from_biopython(bio_feature: SeqFeature, feature: "Feature" = None,
                       leftovers: Dict[str, Any] = None) -> SeqFeature:
        """ Converts a SeqFeature into a single Feature instance.

            Arguments:
                bio_feature: the SeqFeature to convert
                feature: a optional Feature instance to update with the values
                         this class tracks
                leftovers: any qualifiers remaining from the original SeqFeature
                           that have not been used by any subclass

            Returns:
                a Feature instance
        """
        if feature is None:
            feature = Feature(bio_feature.location, bio_feature.type)
            if not leftovers:
                assert isinstance(bio_feature.qualifiers, dict)
                leftovers = bio_feature.qualifiers.copy()
            feature.notes = leftovers.pop("note", [])
        else:
            assert isinstance(feature, Feature)
        if leftovers:
            feature._qualifiers.update(leftovers)
        return feature

    @staticmethod
    def make_qualifiers_copy(bio_feature: SeqFeature) -> Dict[str, Any]:
        """ Makes a shallow copy of a SeqFeature's qualifiers. Only the 'notes'
            key will have a copy taken at a deeper level.
        """
        qualifiers = bio_feature.qualifiers.copy()
        if "note" in qualifiers:
            qualifiers["note"] = qualifiers["note"].copy()
        return qualifiers


class Gene(Feature):
    """ A feature representing a Gene (more general than a CDS) """
    __slots__ = ["_pseudo", "locus_tag", "gene_name"]

    def __init__(self, location, locus_tag=None, gene_name=None, pseudo_gene=False,
                 created_by_antismash=False, qualifiers=None):
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
        return self.locus_tag or self.gene_name

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
    def from_biopython(bio_feature, feature=None, leftovers=None) -> "Gene":
        if leftovers is None:
            leftovers = Feature.make_qualifiers_copy(bio_feature)
        # grab mandatory qualifiers and create the class
        locus = leftovers.pop("locus_tag", [None])[0]
        name = leftovers.pop("gene", [None])[0]
        pseudo = "pseudo" in leftovers
        if pseudo:
            leftovers.pop("pseudo")
        try:
            feature = Gene(bio_feature.location, locus_tag=locus, gene_name=name, pseudo_gene=pseudo)
        except AssertionError:
            print(locus, name, bio_feature.qualifiers)
            raise
        super(Gene, feature).from_biopython(bio_feature, feature=feature, leftovers=leftovers)
        return feature


class ClusterBorder(Feature):
    """ A feature representing a cluster border """
    __slots__ = ["tool", "probability", "cutoff", "extent", "product", "rule",
                 "contig_edge", "high_priority_product"]

    def __init__(self, location: FeatureLocation, tool: str, probability: float = None,
                 cutoff: int = 0, extent: int = 0,
                 product: Optional[str] = None, rule: Optional[str] = None,
                 contig_edge: bool = False, high_priority_product: bool = True) -> None:
        super().__init__(location, feature_type="cluster_border",
                         created_by_antismash=True)
        # required
        self.tool = str(tool)
        # args with simple defaults
        self.high_priority_product = bool(high_priority_product)
        self.contig_edge = bool(contig_edge)
        self.cutoff = int(cutoff)
        self.extent = int(extent)

        # more complicated args
        if product is not None:
            assert isinstance(product, str), type(product)
        self.product = product

        # specific to cluster finder
        self.probability = None
        if probability is not None:
            self.probability = float(probability)

        # specific to rule-based
        if rule is not None:
            assert isinstance(rule, str), type(rule)
        self.rule = rule

    def to_biopython(self, qualifiers=None):
        mine = OrderedDict()
        mine["aStool"] = [self.tool]
        mine["contig_edge"] = [self.contig_edge]
        if self.probability is not None:
            mine["probability"] = [str(self.probability)]
        if self.product:
            mine["product"] = [self.product]
        if self.cutoff:
            mine["cutoff"] = [self.cutoff]
        if self.extent:
            mine["extent"] = [self.extent]
        if self.rule:
            mine["rule"] = [self.rule]
        if qualifiers:
            mine.update(qualifiers)
        return super().to_biopython(mine)

    @staticmethod
    def from_biopython(bio_feature, feature=None, leftovers=None):
        if leftovers is None:
            leftovers = Feature.make_qualifiers_copy(bio_feature)

        # grab mandatory qualifiers and create the class
        tool = leftovers.pop("aStool")[0]

        # optional
        probability = leftovers.pop("probability", [None])[0]
        cutoff = leftovers.pop("cutoff", [0])[0]
        extent = leftovers.pop("extent", [0])[0]
        rule = leftovers.pop("rule", [None])[0]
        product = leftovers.pop("product", [None])[0]
        contig_edge = leftovers.pop("contig_edge", [""])[0] == "True"

        feature = ClusterBorder(bio_feature.location, tool, probability=probability,
                                cutoff=cutoff, extent=extent, rule=rule, product=product,
                                contig_edge=contig_edge)

        # grab parent optional qualifiers
        super(ClusterBorder, feature).from_biopython(bio_feature, feature=feature, leftovers=leftovers)

        return feature

    def __str__(self) -> str:
        return repr(self)

    def __repr__(self) -> str:
        return "ClusterBorder(%s, %s)" % (self.product, self.location)


class AntismashFeature(Feature):
    """ A base class for all sub-CDS Antismash features """
    __slots__ = ["domain_id", "database", "detection", "_evalue", "label",
                 "locus_tag", "_score", "_translation"]

    def __init__(self, location, feature_type):
        super().__init__(location, feature_type, created_by_antismash=True)
        self.domain_id = None
        self.database = None
        self.detection = None
        self._evalue = None  # float
        self.label = None
        self.locus_tag = None
        self._score = None  # float

        self._translation = None

    @property
    def translation(self) -> str:
        """ The amino acid translation of the feature. """
        return self._translation

    @translation.setter
    def translation(self, translation: str):
        self._translation = str(translation)

    @property
    def score(self):
        """ The bitscore reported by a tool when locating the feature """
        return self._score

    @score.setter
    def score(self, score):
        self._score = float(score)

    @property
    def evalue(self):
        """ The e-value reported by a tool when locating the feature """
        return self._evalue

    @evalue.setter
    def evalue(self, evalue):
        self._evalue = float(evalue)

    def to_biopython(self, qualifiers=None):
        mine = OrderedDict()
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
    def from_biopython(bio_feature, qualifiers=None, feature=None, leftovers=None):
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
        feature.translation = leftovers.pop("translation", [None])[0]
        if "evalue" in leftovers:
            feature.evalue = float(leftovers.pop("evalue")[0])
        if "score" in leftovers:
            feature.score = float(leftovers.pop("score")[0])

        # grab parent optional qualifiers
        return Feature.from_biopython(bio_feature, feature=feature, leftovers=leftovers)


class Domain(AntismashFeature):
    """ A base class for features which represent a domain type """
    __slots__ = ["tool", "domain"]

    def __init__(self, location, feature_type):
        super().__init__(location, feature_type)
        self.tool = None
        self.domain = None

    def to_biopython(self, qualifiers=None):
        mine = OrderedDict()
        if self.tool:
            mine["aSTool"] = [self.tool]
        if self.domain:
            mine["aSDomain"] = [self.domain]
        if self.domain_id:
            mine["aSDomain_id"] = [self.domain_id]
        if qualifiers:
            mine.update(qualifiers)
        return super().to_biopython(mine)

    @staticmethod
    def from_biopython(bio_feature, feature=None, leftovers=None):
        if leftovers is None:
            leftovers = Feature.make_qualifiers_copy(bio_feature)
        if not feature:
            raise ValueError("Domain shouldn't be instantiated directly")
        else:
            assert isinstance(feature, Domain), type(feature)

        # grab optional qualifiers
        if "aSTool" in leftovers:
            feature.tool = leftovers.pop("aSTool")[0]
        if "asDomain" in leftovers:
            feature.domain = leftovers.pop("asDomain")[0]

        # grab parent optional qualifiers
        return AntismashFeature.from_biopython(bio_feature, feature=feature, leftovers=leftovers)


class CDSMotif(Domain):
    """ A base class for features that represent a motif within a CDSFeature """
    __slots__ = ["motif"]

    def __init__(self, location):
        super().__init__(location, feature_type="CDS_motif")
        self.motif = None

    @staticmethod
    def from_biopython(bio_feature, feature=None, leftovers=None):
        if leftovers is None:
            leftovers = Feature.make_qualifiers_copy(bio_feature)
        if not feature:
            feature = CDSMotif(bio_feature.location)

        feature.motif = leftovers.pop("description", [None])[0]
        return super(CDSMotif, feature).from_biopython(bio_feature, feature, leftovers)

    def to_biopython(self, qualifiers=None):
        mine = OrderedDict()
        if self.motif:
            mine["motif"] = [self.motif]
        if qualifiers:
            mine.update(qualifiers)
        return super().to_biopython(mine)


class PFAMDomain(Domain):
    """ A feature representing a PFAM domain within a CDS.
    """
    __slots__ = ["description", "db_xref", "probability"]

    def __init__(self, location: FeatureLocation, description: str, domain: Optional[str] = None) -> None:
        super().__init__(location, feature_type="PFAM_domain")
        assert description and isinstance(description, str), description
        if domain is not None:
            assert domain and isinstance(domain, str), domain
        self.domain = domain
        self.description = description
        self.probability = None
        self.db_xref = []  # type: List[str]

    def to_biopython(self, qualifiers=None):
        mine = OrderedDict()
        mine["description"] = self.description
        if self.probability is not None:
            mine["probability"] = [self.probability]
        if self.db_xref:
            mine["db_xref"] = self.db_xref
        if self.domain is not None:
            mine["domain"] = self.domain
        if qualifiers:
            mine.update(qualifiers)
        return super().to_biopython(mine)

    @staticmethod
    def from_biopython(bio_feature, feature=None, leftovers=None):
        if leftovers is None:
            leftovers = Feature.make_qualifiers_copy(bio_feature)
        # grab mandatory qualifiers and create the class
        description = leftovers.pop("description")[0]
        feature = PFAMDomain(bio_feature.location, description)

        # grab optional qualifiers
        feature.db_xref = leftovers.pop("db_xref", [])

        # grab parent optional qualifiers
        return super(PFAMDomain, feature).from_biopython(bio_feature, feature=feature, leftovers=leftovers)


class AntismashDomain(Domain):
    """ A class to represent a Domain with extra specificities and type information """
    __slots__ = ["domain_subtype", "specificity"]

    def __init__(self, location):
        super().__init__(location, feature_type="aSDomain")
        self.domain_subtype = None
        self.specificity = []

    def to_biopython(self, qualifiers=None) -> List[SeqFeature]:
        mine = OrderedDict()  # type: Dict[str, List[str]]
        if self.domain_subtype:
            mine["domain_subtype"] = [self.domain_subtype]
        if self.specificity:
            mine["specificity"] = self.specificity
        if qualifiers:
            mine.update(qualifiers)
        return super().to_biopython(mine)

    @staticmethod
    def from_biopython(bio_feature, feature=None, leftovers=None) -> "AntismashDomain":
        if leftovers is None:
            leftovers = Feature.make_qualifiers_copy(bio_feature)
        # grab mandatory qualifiers and create the class
        feature = AntismashDomain(bio_feature.location)

        # grab optional qualifiers
        feature.domain_subtype = leftovers.pop("domain_subtype", [None])[0]
        feature.specificity = list(leftovers.pop("specificity", []))

        # grab parent optional qualifiers
        super(AntismashDomain, feature).from_biopython(bio_feature, feature=feature, leftovers=leftovers)

        return feature


class CDSFeature(Feature):
    """ A feature representing a single CDS/gene. """
    __slots__ = ["_translation", "protein_id", "locus_tag", "gene", "product",
                 "transl_table", "_sec_met", "product_prediction", "cluster", "_gene_functions",
                 "unique_id", "_nrps_pks", "motifs"]

    def __init__(self, location, translation=None, locus_tag=None, protein_id=None,
                 product=None, gene=None):
        super().__init__(location, feature_type="CDS")
        # mandatory
        #  codon_start
        #  db_xref
        self._gene_functions = GeneFunctionAnnotations()

        # semi-optional
        self.protein_id = _sanitise_id_value(protein_id)
        self.locus_tag = _sanitise_id_value(locus_tag)
        self.gene = _sanitise_id_value(gene)
        self._translation = None
        if translation is not None:
            self.translation = translation

        # optional
        self.product = product
        self.transl_table = "Standard"
        self._sec_met = None  # SecMetQualifier()
        self._nrps_pks = NRPSPKSQualifier()
        self.product_prediction = []  # TODO: shift into nrps sub section?

        self.motifs = []

        if not (protein_id or locus_tag or gene):
            raise ValueError("CDSFeature requires at least one of: gene, protein_id, locus_tag")

        # runtime-only data
        self.cluster = None
        self.unique_id = None  # set only when added to a record

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
    def sec_met(self, sec_met):
        if sec_met is not None and not isinstance(sec_met, SecMetQualifier):
            raise TypeError("CDSFeature.sec_met can only be set to an instance of SecMetQualifier")
        self._sec_met = sec_met

    @property
    def nrps_pks(self):
        """ The NRPSPKSQualifier of the feature """
        return self._nrps_pks

    @nrps_pks.setter
    def nrps_pks(self, qualifier):
        if qualifier is not None and not isinstance(qualifier, NRPSPKSQualifier):
            raise TypeError("CDSFeature.nrps_pks can only be set to an instance of NRPSPKSQualifier")
        self._nrps_pks = qualifier

    @property
    def translation(self) -> str:
        """ The translation of the CDS, as a string of amino acids """
        return self._translation

    @translation.setter
    def translation(self, translation: str) -> None:
        assert "-" not in translation, "%s contains - in translation" % self.get_name()
        self._translation = str(translation)

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
    def from_biopython(bio_feature: SeqFeature, feature: Feature = None,
                       leftovers: Optional[Dict] = None) -> "CDSFeature":
        if leftovers is None:
            leftovers = Feature.make_qualifiers_copy(bio_feature)
        # grab mandatory qualifiers and create the class

        # semi-optional qualifiers
        protein_id = leftovers.pop("protein_id", [None])[0]
        locus_tag = leftovers.pop("locus_tag", [None])[0]
        gene = leftovers.pop("gene", [None])[0]
        if not (gene or protein_id or locus_tag):
            if "pseudo" in leftovers or "pseudogene" in leftovers:
                locus_tag = "pseudo_%d" % int(bio_feature.location.start + 1)  # 1-indexed
            else:
                # TODO solve somehow?
                logging.critical("CDS feature created from biopython without identifier: %s", bio_feature)
                raise ValueError("CDSFeature requires at least one of: gene, protein_id, locus_tag")
        translation = leftovers.pop("translation", [None])[0]
        if translation and "-" in translation:
            logging.warning("Translation for CDS %s (at %s) has a gap. Discarding and regenerating.",
                            locus_tag or protein_id or gene, bio_feature.location)
            translation = None

        feature = CDSFeature(bio_feature.location, translation, gene=gene,
                             locus_tag=locus_tag, protein_id=protein_id)

        # grab optional qualifiers
        feature.product = leftovers.pop("product", [None])[0]
        feature.transl_table = leftovers.pop("transl_table", ["Standard"])[0]
        sec_met = leftovers.pop("sec_met", None)
        if sec_met:
            feature.sec_met = SecMetQualifier.from_biopython(sec_met)
        gene_functions = leftovers.pop("gene_functions", [])
        if gene_functions:
            feature.gene_functions.add_from_qualifier(gene_functions)
        feature.product_prediction = leftovers.pop("aSProdPred", [])

        # grab parent optional qualifiers
        super(CDSFeature, feature).from_biopython(bio_feature, feature=feature, leftovers=leftovers)

        return feature

    def to_biopython(self, qualifiers: Dict[str, List[str]] = None) -> SeqFeature:
        mine = OrderedDict()  # type: Dict[str, List[str]]
        # mandatory
        mine["translation"] = [self.translation]
        if self.product_prediction:
            mine["aSProdPred"] = [self.product_prediction]
        # optional
        for attr in ["gene", "transl_table", "locus_tag",
                     "protein_id", "product"]:
            val = getattr(self, attr)
            if val:
                mine[attr] = [str(val)]
        if self._gene_functions:
            mine["gene_functions"] = list(map(str, self._gene_functions))
        # since it's already a list
        if self.sec_met:
            mine["sec_met"] = self.sec_met
        # respect qualifiers given to us
        if qualifiers:
            mine.update(qualifiers)
        return super().to_biopython(mine)

    def __repr__(self) -> str:
        return str(self)

    def __str__(self) -> str:
        return "CDS(%s, %s)" % (self.get_name(), self.location)


class Prepeptide(CDSMotif):
    """ A class representing a prepeptide. Used for tracking a multi-feature
        construction with a leader, core and tail. To allow for multiple types
        of prepeptide (e.g. lanthi- or sacti-peptides), only the core must exist.
    """
    def __init__(self, location, peptide_class, core, locus_tag, peptide_subclass=None,
                 score=0., monoisotopic_mass=0., molecular_weight=0.,
                 alternative_weights=None, leader="", tail="", **kwargs):
        """
            Arguments:
                peptide_class: the kind of prepeptide, e.g. 'lanthipeptide', 'thiopeptide'
                core: the sequence of the core
                locus_tag: the locus tag to use for the feature
                prepeptide_subclass: the subclass of the prepeptide, e.g. 'Type II'
                leader: the sequence of the leader, if it exists
                tail: the sequence of the tail, if it exists
                ... other args that will be passed through to CDSFeature
        """
        for arg in [peptide_class, core, leader, tail]:
            assert isinstance(arg, str), type(arg)
        self._leader = leader
        self._core = core
        self._tail = tail
        super().__init__(location, **kwargs)
        self.locus_tag = locus_tag
        self.type = "CDS_motif"
        self.peptide_class = peptide_class
        if peptide_subclass:
            peptide_subclass = peptide_subclass.replace("-", " ")  # "Type-II" > "Type II"
        self.peptide_subclass = peptide_subclass
        self.score = float(score)
        self.monoisotopic_mass = float(monoisotopic_mass)
        self.molecular_weight = float(molecular_weight)
        self.alternative_weights = []
        if alternative_weights:
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


class Cluster(Feature):
    """ A feature representing a cluster. Tracks which CDS features belong to it"""
    __slots__ = ["_extent", "_cutoff", "_products", "contig_edge",
                 "detection_rules", "smiles_structure",
                 "clusterblast", "knownclusterblast", "subclusterblast",
                 "parent_record", "cds_children", "borders", "monomers_prediction"]

    def __init__(self, location: FeatureLocation, cutoff: int, extent: int, products: List) -> None:
        super().__init__(location, feature_type="cluster",
                         created_by_antismash=True)

        self._extent = int(extent)
        self._cutoff = int(cutoff)
        self._products = []
        for product in products:
            self.add_product(product)

        self.contig_edge = None  # hmm_detection borderpredict
        self.detection_rules = []
        self.smiles_structure = None  # SMILES string
        self.monomers_prediction = None

        self.clusterblast = None
        self.knownclusterblast = None
        self.subclusterblast = None

        # for runtime management
        self.parent_record = None
        self.cds_children = OrderedDict()
        self.borders = []

    @property
    def products(self) -> Iterable[str]:
        """ The products of a cluster """
        return tuple(self._products)

    def add_product(self, product: str) -> None:
        """ Add the given product to the cluster's list of products """
        assert product and isinstance(product, str), str(product)
        self._products.append(product)

    def get_cluster_number(self):
        """ Returns the cluster number which the parent record uses to refer to
            this cluster. """
        if not self.parent_record:
            raise ValueError("Cluster not contained in record")
        return self.parent_record.get_cluster_number(self)

    def trim_overlapping(self):
        """ Shrinks the cluster, where possible, to exclude any features which
            overlap with the edges of the cluster.
            Any feature fully contained before shrinking will still be fully
            contained.
        """
        if not self.parent_record:
            logging.warning("Trimming cluster which does not belong to a record")
            return
        features = self.parent_record.get_cds_features_within_location(self.location,
                                            with_overlapping=True)
        # don't trim if there's no features to trim by
        if not features:
            return

        # find the deepest feature that only overlaps at the beginning
        previous = None
        index = 0
        current = features[index]
        # track where to trim to
        start = self.location.start
        while current.overlaps_with(self) and not current.is_contained_by(self):
            start = max([start, current.location.start, current.location.end])
            previous = current
            index += 1
            if index >= len(features):
                current = None
                break
            current = features[index]

        # don't cause a contained feature to now overlap only
        if previous and current:
            start = min([start, current.location.start, current.location.end])

        # find the deepest feature that only overlaps at the end
        # but skip any indices already covered in the lead search
        lead_index = index
        previous = None
        index = len(features) - 1
        current = features[index]
        # track where to trim to
        end = self.location.end
        while index > lead_index and current.overlaps_with(self) and not current.is_contained_by(self):
            end = min([end, current.location.start, current.location.end])
            previous = current
            index -= 1
            if index < 0:
                current = None
                break
            current = features[index]

        # but don't cause a contained feature to now overlap only
        if previous and current:
            end = max([end, current.location.start, current.location.end])

        # finally, do the trim itself
        new_loc = FeatureLocation(start, end, self.location.strand)
        if self.location.start != start or self.location.end != end:
            logging.debug("Cluster %d trimming location from %s to %s",
                          self.get_cluster_number(), self.location, new_loc)
        # make sure the size is never increased
        assert self.location.start <= start < end <= self.location.end
        self.location = new_loc

        for cds in self.cds_children:
            assert cds.is_contained_by(self), "cluster trimming removed wholly contained CDS"

    def add_cds(self, cds: CDSFeature):
        """ Adds a CDSFeature to the cluster """
        assert isinstance(cds, CDSFeature)
        assert cds.is_contained_by(self), "cds %s outside cluster %s" % (cds, self)
        self.cds_children[cds] = None

    @property
    def cutoff(self):
        """ The maximal distance between genes when defining the cluster.
            The distance between core genes after definition will likely be
            smaller."""
        return self._cutoff

    @cutoff.setter
    def cutoff(self, cutoff):
        if cutoff is not None:
            cutoff = int(cutoff)
        self._cutoff = cutoff

    @property
    def extent(self):
        """ The distance the cluster extends from the first and last genes which
            from which the cluster was defined.
        """
        return self._extent

    @extent.setter
    def extent(self, extent):
        if extent is not None:
            self._extent = int(extent)

    def get_product_string(self):
        """ Returns the cluster's products as a single string """
        assert None not in self._products, self._products
        return "-".join(self._products)

    @property
    def probability(self) -> Optional[float]:
        """ The cluster probability, if relevant. """
        probabilities = {border.probability for border in self.borders}
        # one border ignores probabilities, then don't use a probability
        if None in probabilities:
            return None
        # if all agree on the probability
        if len(probabilities) == 1:
            return list(probabilities)[0]
        # if they disagree, return None again
        return None

    @staticmethod
    def from_biopython(bio_feature, feature=None, leftovers=None):
        if leftovers is None:
            leftovers = Feature.make_qualifiers_copy(bio_feature)
        # grab mandatory qualifiers and create the class
        cutoff = leftovers.pop("cutoff")[0]
        extent = leftovers.pop("extension")[0]
        products = leftovers.pop("product")[0].split("-")
        cluster = Cluster(bio_feature.location, cutoff, extent, products)
        # take the detection rules from "note"
        # first check it exists
        index = -1
        for i, note in enumerate(leftovers.get("note", [])):
            if note.startswith("Detection rule"):
                index = i
                break
        # then split it into the relevant pieces
        if index > -1:
            products = []
            rules = []
            text = leftovers["note"].pop(index)
            text = text.split(":", 1)[1]  # strip the leadin
            text = text.split(";")  # separate rules
            for rule in text:
                if not rule:
                    continue
                assert ": " in rule, rule
                product, rule = rule.split(": ")
                rule = rule[1:-1]  # strip ( )
                products.append(product.strip())
                rules.append(rule)
            cluster.detection_rules = rules

        # grab optional qualifiers
        contig_edge = leftovers.pop("contig_edge", [None])[0]
        if not contig_edge:
            cluster.contig_edge = None
        else:
            cluster.contig_edge = contig_edge == "True"
        cluster.smiles_structure = leftovers.pop("structure", None)
        # grab optional parent qualifiers
        super(Cluster, cluster).from_biopython(bio_feature, cluster, leftovers)

        return cluster

    def to_biopython(self, qualifiers=None):
        mine = OrderedDict()
        mine["cutoff"] = [str(self.cutoff)]
        mine["extension"] = [str(self.extent)]
        if self.contig_edge is not None:
            mine["contig_edge"] = [str(self.contig_edge)]
        assert isinstance(self._products, list), type(self.products)
        mine["product"] = [self.get_product_string()]
        rule_text = ["Detection rule(s) for this cluster type:"]
        assert isinstance(self.detection_rules, list), type(self.detection_rules)
        for product, rule in zip(self.products, self.detection_rules):
            rule_text.append("%s: (%s);" % (product, rule))
        rule_text = " ".join(rule_text)
        if qualifiers:
            mine.update(qualifiers)
        if "note" not in mine:
            mine["note"] = []
        mine["note"].append(rule_text)
        return super().to_biopython(mine)

    def write_to_genbank(self, filename=None, directory=None, record=None):
        """ Writes a genbank file containing only the information contained
            within the Cluster.
        """
        if not filename:
            filename = "%s.cluster%03d.gbk" % (self.parent_record.id, self.get_cluster_number())
        if directory:
            filename = os.path.join(directory, filename)

        if record is None:
            record = self.parent_record.to_biopython()
        assert isinstance(record, SeqRecord)
        with warnings.catch_warnings():
            warnings.simplefilter("ignore")
            cluster_record = record[self.location.start:self.location.end]

        cluster_record.annotations["date"] = record.annotations.get("date", '')
        cluster_record.annotations["source"] = record.annotations.get("source", '')
        cluster_record.annotations["organism"] = record.annotations.get("organism", '')
        cluster_record.annotations["taxonomy"] = record.annotations.get("taxonomy", [])
        cluster_record.annotations["data_file_division"] = record.annotations.get("data_file_division", 'UNK')
        # our cut-out clusters are always linear
        cluster_record.annotations["topology"] = "linear"

        seqio.write([cluster_record], filename, 'genbank')


def _sanitise_id_value(name: Optional[str]) -> Optional[str]:
    """ Ensures a name doesn't contain characters that will break external programs"""
    if name is None:
        return None
    name = str(name)
    illegal_chars = set("!\"#$%&()*+,:; \r\n\t=>?@[]^`'{|}/ ")
    for char in set(name).intersection(illegal_chars):
        name = name.replace(char, "_")
    return name
