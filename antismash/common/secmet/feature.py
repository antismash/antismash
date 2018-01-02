# License: GNU Affero General Public License v3 or later
# A copy of GNU AGPL v3 should have been included in this software package in LICENSE.txt.

from collections import defaultdict, OrderedDict
from enum import Enum, unique
import logging
import os
import warnings
from typing import Any, Dict, List, Optional, Tuple, Union

from helperlibs.bio import seqio

from Bio.Seq import Seq
from Bio.SeqFeature import SeqFeature, FeatureLocation, CompoundLocation
from Bio.SeqRecord import SeqRecord

from .qualifiers import NRPSPKSQualifier, SecMetQualifier


class Feature:
    """ The base class of any feature. Contains only a location, the label of the
        subclass, the 'notes' qualifier, and other qualifiers not tracked by any
        subclass.
    """
    __slots__ = ["location", "notes", "type", "_qualifiers"]

    def __init__(self, location: FeatureLocation, feature_type: str) -> None:
        assert isinstance(location, (FeatureLocation, CompoundLocation)), type(location)
        self.location = location
        self.notes = []  # type: List[str]
        assert feature_type
        self.type = str(feature_type)
        self._qualifiers = OrderedDict()  # type: Dict[str, List[str]]

    @property
    def strand(self) -> int:
        """ A simple wrapper to access the location strand directly.
        """
        return self.location.strand

    def extract(self, sequence: Union[Seq, str]) -> Union[Seq, str]:
        """ Extracts the section of the given sequence that this feature covers.

            Return type is the same as the input type.
        """
        return self.location.extract(sequence)

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
        return self.location.start in location \
               or self.location.end - 1 in location \
               or location.start in self.location \
               or location.end - 1 in self.location

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


class ClusterBorder(Feature):
    """ A feature representing a cluster border """
    __slots__ = ["tool", "probability", "parent_cluster"]

    def __init__(self, location, tool, probability):
        super().__init__(location, feature_type="cluster_border")
        self.notes.append("best prediction")
        self.tool = str(tool)
        self.probability = float(probability)

        # runtime only
        self.parent_cluster = None

    def to_biopython(self, qualifiers=None):
        mine = OrderedDict()
        mine["tool"] = [self.tool]
        mine["probability"] = [str(self.probability)]
        if qualifiers:
            mine.update(qualifiers)
        return super().to_biopython(mine)

    @staticmethod
    def from_biopython(bio_feature, feature=None, leftovers=None):
        if leftovers is None:
            leftovers = Feature.make_qualifiers_copy(bio_feature)

        # grab mandatory qualifiers and create the class
        tool = leftovers.pop("tool")[0]
        # clean up the note so it isn't duplicated
        notes = leftovers.get("note")
        if notes:
            for i, note in enumerate(notes):
                if note == "best prediction":
                    notes.pop(i)
                    break
        probability = leftovers.pop("probability")[0]
        feature = ClusterBorder(bio_feature.location, tool, probability)

        # grab parent optional qualifiers
        super(ClusterBorder, feature).from_biopython(bio_feature, feature=feature, leftovers=leftovers)

        return feature


class AntismashFeature(Feature):
    """ A base class for all sub-CDS Antismash features """
    __slots__ = ["domain_id", "database", "detection", "_evalue", "label",
                 "locus_tag", "_score", "_translation"]

    def __init__(self, location, feature_type):
        super().__init__(location, feature_type)
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


class Domain(AntismashFeature):
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


class CDSMotif(Domain):
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
    __slots__ = ["description", "db_xref"]

    def __init__(self, location: FeatureLocation, description: str) -> None:
        super().__init__(location, feature_type="PFAM_domain")
        assert isinstance(description, str)
        self.description = description
        self.db_xref = []  # type: List[str]

    def to_biopython(self, qualifiers=None):
        mine = OrderedDict()
        mine["description"] = self.description
        if self.db_xref:
            mine["db_xref"] = self.db_xref
        if qualifiers:
            mine.update(qualifiers)
        return super().to_biopython(mine)

    @staticmethod
    def from_biopython(bio_feature, feature=None, leftovers=None):
        if leftovers is None:
            leftovers = Feature.make_qualifiers_copy(bio_feature)
        # grab mandatory qualifiers and create the class
        description = leftovers.pop("description")
        feature = PFAMDomain(bio_feature.location, description)

        # grab optional qualifiers
        feature.db_xref = leftovers.pop("db_xref", [])

        # grab parent optional qualifiers
        super(PFAMDomain, feature).from_biopython(bio_feature, feature=feature, leftovers=leftovers)

        return feature


class AntismashDomain(Domain):
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


@unique
class GeneFunction(Enum):
    """ An Enum representing the function of a gene.
        Allows for more flexible conversion and more robust value constraints.
    """
    OTHER = 0
    CORE = 1
    ADDITIONAL = 2
    TRANSPORT = 3
    REGULATORY = 4

    def __str__(self) -> str:
        # because this information ends up in the record, make these more
        # labels more meaningful for users
        if self == GeneFunction.CORE:
            return "biosynthetic"
        if self == GeneFunction.ADDITIONAL:
            return "biosynthetic-additional"
        return str(self.name).lower()

    @staticmethod
    def from_string(label: str) -> "GeneFunction":
        """ Converts a string to a GeneFunction instance when possible.
            Raises an error if not possible.
        """
        for value in GeneFunction:
            if str(value) == label:
                return value
        raise ValueError("Unknown gene function label: %s", label)


class GeneFunctionAnnotations:
    """ Tracks multiple annotations of gene functions. Each annotation contains
        the declared function of the gene, which tool determined the function,
        and a comment on how the determination was made.
    """
    slots = ["_annotations", "_by_tool", "_by_function"]

    class GeneFunctionAnnotation:
        """ A single instance of an annotation. """
        slots = ["function", "tool", "description"]

        def __init__(self, function: GeneFunction, tool: str, description: str) -> None:
            assert isinstance(function, GeneFunction), "wrong type: %s" % type(function)
            assert tool and len(tool.split()) == 1, tool  # no whitespace in tool name
            assert description
            self.function = function
            self.tool = tool
            self.description = description

        def __str__(self):
            return "%s (%s) %s" % (self.function, self.tool, self.description)

        def __repr__(self):
            return "GeneFunctionAnnotation(function=%r, tool='%s', '%s')" % (self.function, self.tool, self.description)

    def __init__(self):
        self._annotations = []
        self._by_tool = defaultdict(list)
        self._by_function = defaultdict(list)

    def __iter__(self):
        for annotation in self._annotations:
            yield annotation

    def __len__(self):
        return len(self._annotations)

    def add(self, function: GeneFunction, tool: str, description: str) -> "GeneFunctionAnnotations.GeneFunctionAnnotation":
        """ Adds a gene function annotation. If an existing annotation has all
            the same values as those provided, a duplicate will not be added.

            Arguments:
                function: a GeneFunction value
                tool: a string naming the tool that determined the gene function
                description: description text for storing extra information

            Returns:
                the GeneFunction added (or the existing one if duplicated)
        """
        # if there's already an exactly similar function annotation, skip adding
        existing_functions = self._by_function.get(function, [])
        for existing in existing_functions:
            if existing.tool == tool and existing.description == description:
                return existing
        new = GeneFunctionAnnotations.GeneFunctionAnnotation(function, tool, description)
        self._by_tool[tool].append(new)
        self._by_function[function].append(new)
        self._annotations.append(new)
        return new

    def add_from_qualifier(self, qualifier: List[str]):
        """ Converts a string-based qualifier into one or more GeneFunctions and
            adds them to the current set.
            Expected format will be as per str(GeneFunctionAnnotation).
        """
        for section in qualifier:
            function, tool, description = section.split(maxsplit=2)
            tool = tool[1:-1]  # strip the ()
            self.add(GeneFunction.from_string(function), tool, description)

    def get_by_tool(self, tool: str) -> Optional[List]:
        """ Returns a list of all GeneFunctionAnnotations which were added with
            the tool name provided.
        """
        return self._by_tool.get(tool)

    def get_by_function(self, function: GeneFunction) -> Optional[List]:
        """ Returns a list of GeneFunctionAnnotations which have the same
            gene function as that provided.
        """
        return self._by_function.get(function)

    def get_classification(self) -> GeneFunction:
        """ Returns the function of the gene, in the priority order of priority:
             - CORE if cluster defined by this gene
             - the function determined by smCOGs, if it exists
             - a function from all tools if they all agree
             - or OTHER
        """
        # if no annotations, skip to OTHER
        if not self._annotations:
            return GeneFunction.OTHER
        # priority for tools
        for tool in ["cluster_definition", "smcogs"]:
            annotations = self._by_tool[tool]
            if annotations:
                return annotations[0].function
        # otherwise check all agree
        function = self._annotations[0].function
        for annotation in self._annotations[1:]:
            if annotation.function != function:
                return GeneFunction.OTHER
        return function


class CDSFeature(Feature):
    """ A feature representing a single CDS/gene. """
    __slots__ = ["_translation", "protein_id", "locus_tag", "gene", "product",
                 "transl_table", "_sec_met", "aSProdPred", "cluster", "_gene_functions",
                 "unique_id", "_nrps_pks", "motifs"]

    def __init__(self, location, translation=None, locus_tag=None, protein_id=None,
                 product=None, gene=None):
        super().__init__(location, feature_type="CDS")
        # mandatory
        #  codon_start
        #  db_xref
        self._gene_functions = GeneFunctionAnnotations()

        # semi-optional
        self._translation = None
        if translation is not None:
            self.translation = translation
        self.protein_id = protein_id
        self.locus_tag = locus_tag
        self.gene = gene

        # optional
        self.product = product
        self.transl_table = None
        self._sec_met = None  # SecMetQualifier()
        self._nrps_pks = NRPSPKSQualifier()
        self.aSProdPred = []  # TODO: shift into nrps sub section?

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
        assert "-" not in translation
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
        translation = leftovers.pop("translation", [None])[0]
        protein_id = leftovers.pop("protein_id", [None])[0]
        locus_tag = leftovers.pop("locus_tag", [None])[0]
        gene = leftovers.pop("gene", [None])[0]
        if not (gene or protein_id or locus_tag):
            # TODO solve somehow
            logging.critical("CDS feature created from biopython without identifier")
            raise ValueError("CDSFeature requires at least one of: gene, protein_id, locus_tag")

        feature = CDSFeature(bio_feature.location, translation, gene=gene,
                             locus_tag=locus_tag, protein_id=protein_id)

        # grab optional qualifiers
        feature.product = leftovers.pop("product", [None])[0]
        feature.transl_table = leftovers.pop("transl_table", [None])[0]
        sec_met = leftovers.pop("sec_met", None)
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
        for attr in ["gene", "aSProdPred", "transl_table", "locus_tag",
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
        self.peptide_subclass = peptide_subclass.replace("-", " ")  # "Type-II" > "Type II"
        self.score = float(score)
        self.monoisotopic_mass = float(monoisotopic_mass)
        self.molecular_weight = float(molecular_weight)
        if alternative_weights:
            self.alternative_weights = [float(weight) for weight in alternative_weights]

    @property
    def translation(self) -> str:
        return self._leader + self._core + self._tail

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
    __slots__ = ["_extent", "_cutoff", "products", "contig_edge",
                 "detection_rules", "smiles_structure", "probability",
                 "clusterblast", "knownclusterblast", "subclusterblast",
                 "parent_record", "cds_children", "borders", "monomers_prediction"]

    def __init__(self, location, cutoff, extent, products):
        super().__init__(location, feature_type="cluster")

        self._extent = int(extent)
        self._cutoff = int(cutoff)
        assert isinstance(products, list)
        self.products = products

        self.contig_edge = None  # hmm_detection borderpredict
        self.detection_rules = []
        self.smiles_structure = None  # SMILES string
        self.monomers_prediction = None
        self.probability = None  # clusterfinder probability # TODO: unify with notes version

        self.clusterblast = None
        self.knownclusterblast = None
        self.subclusterblast = None

        # for runtime management
        self.parent_record = None
        self.cds_children = OrderedDict()
        self.borders = []

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
        logging.debug("Cluster %d trimming location from %s to %s", self.get_cluster_number(),
                      self.location, new_loc)
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
        return "-".join(self.products)

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
            assert sorted(products) == sorted(cluster.products)
            cluster.products = products
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
        assert isinstance(self.products, list), type(self.products)
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
            record = self.parent_record.to_biopython
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
