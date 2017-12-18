# License: GNU Affero General Public License v3 or later
# A copy of GNU AGPL v3 should have been included in this software package in LICENSE.txt.

from collections import defaultdict, OrderedDict
from enum import Enum, unique
import logging
import os
import warnings
from typing import Dict, List, Optional, Tuple, Union

from helperlibs.bio import seqio

from Bio.SeqFeature import SeqFeature, FeatureLocation, CompoundLocation
from Bio.SeqRecord import SeqRecord


class Feature:
    __slots__ = ["location", "notes", "type", "_qualifiers"]

    def __init__(self, location, feature_type):
        assert isinstance(location, (FeatureLocation, CompoundLocation)), type(location)
        self.location = location
        self.notes = []
        assert feature_type
        self.type = str(feature_type)
        self._qualifiers = OrderedDict()

    @property
    def strand(self):
        return self.location.strand

    def extract(self, sequence):
        return self.location.extract(sequence)

    def get_qualifier(self, key):
        qualifier = self._qualifiers.get(key)
        if qualifier:
            return tuple(qualifier)
        return None

    def overlaps_with(self, other: Union["Feature", FeatureLocation]):
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

    def is_contained_by(self, other: Union["Feature", FeatureLocation]):
        end = self.location.end - 1  # to account for the non-inclusive end
        if isinstance(other, Feature):
            return self.location.start in other.location and end in other.location
        if isinstance(other, FeatureLocation):
            return self.location.start in other and end in other
        raise TypeError("Container must be a Feature or a FeatureLocation, not %s" % type(other))

    def to_biopython(self, qualifiers=None):
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

    def __lt__(self, other):
        """ Allows sorting Features by location without key complication """
        assert isinstance(other, Feature)
        if self.location.start < other.location.start:
            return True
        elif self.location.start == other.location.start:
            return self.location.end < other.location.end
        return False

    def __str__(self):
        return repr(self)

    def __repr__(self):
        return "%s(%s)" % (self.type, self.location)

    @staticmethod
    def from_biopython(bio_feature, feature=None, leftovers=None):
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
    def make_qualifiers_copy(bio_feature):
        qualifiers = bio_feature.qualifiers.copy()
        if "note" in qualifiers:
            qualifiers["note"] = qualifiers["note"].copy()
        return qualifiers


class ClusterBorder(Feature):
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
    __slots__ = ["domain_id", "database", "detection", "_evalue", "label",
                 "locus_tag", "_score", "translation"]

    def __init__(self, location, feature_type):
        super().__init__(location, feature_type)
        self.domain_id = None
        self.database = None
        self.detection = None
        self._evalue = None  # float
        self.label = None
        self.locus_tag = None
        self._score = None  # float

        self.translation = None

    @property
    def score(self):
        return self._score

    @score.setter
    def score(self, score):
        self._score = float(score)

    @property
    def evalue(self):
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
        if self.translation:
            mine["translation"] = [self.translation]
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
    __slots__ = ["tool", "domain", "domain_id"]

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
    __slots__ = ["description", "db_xref"]

    def __init__(self, location: FeatureLocation, description: str) -> None:
        super().__init__(location, feature_type="PFAM_domain")
        assert isinstance(description, str)
        self.description = description
        self.db_xref = []

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

    def to_biopython(self, qualifiers=None):
        mine = OrderedDict()
        if self.domain_subtype:
            mine["domain_subtype"] = [self.domain_subtype]
        if self.specificity:
            mine["specificity"] = self.specificity
        if qualifiers:
            mine.update(qualifiers)
        return super().to_biopython(mine)

    @staticmethod
    def from_biopython(bio_feature, feature=None, leftovers=None):
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
    OTHER = 0
    CORE = 1
    ADDITIONAL = 2
    TRANSPORT = 3
    REGULATORY = 4

    def __str__(self):
        # because this information ends up in the record, make these more
        # labels more meaningful for users
        if self == GeneFunction.CORE:
            return "biosynthetic"
        if self == GeneFunction.ADDITIONAL:
            return "biosynthetic-additional"
        return str(self.name).lower()

    @staticmethod
    def from_string(label):
        for value in GeneFunction:
            if str(value) == label:
                return value
        raise ValueError("Unknown gene function label: %s", label)


class GeneFunctionAnnotations:
    slots = ["_annotations", "_by_tool", "_by_function"]

    class GeneFunctionAnnotation:
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

    def add_from_qualifier(self, qualifier: List):
        for section in qualifier:
            function, tool, description = section.split(maxsplit=2)
            tool = tool[1:-1]  # strip the ()
            self.add(GeneFunction.from_string(function), tool, description)

    def get_by_tool(self, tool: str) -> Optional[List]:
        return self._by_tool.get(tool)

    def get_by_function(self, function: GeneFunction) -> Optional[List]:
        return self._by_function.get(function)

    def get_classification(self) -> GeneFunction:
        """ Returns, in order of priority:
             - CORE if cluster defined by this gene
             - smCOGs function if it exists
             - a function from other tools if they all agree
             - OTHER
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

        # TODO: add when active site finder completed
        #self.aSASF_choice = []
        #self.aSASF_note = []
        #self.aSASF_prediction = []
        #self.aSASF_scaffold = []

        # runtime-only data
        self.cluster = None
        self.unique_id = None  # set only when added to a record

    @property
    def gene_functions(self):
        return self._gene_functions

    @property
    def gene_function(self):
        return self._gene_functions.get_classification()

    @property
    def sec_met(self):
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
        return self._translation

    @translation.setter
    def translation(self, translation) -> None:
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


class Prepeptide(CDSFeature):
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
        super().__init__(location, locus_tag=locus_tag, **kwargs)
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
        return self._leader

    @leader.setter
    def leader(self, leader: str) -> None:
        assert isinstance(leader, str)
        self._leader = leader

    @property
    def core(self) -> str:
        return self._core

    @leader.setter
    def leader(self, core: str) -> None:
        assert isinstance(core, str)
        self._core = core

    @property
    def tail(self) -> str:
        return self._tail

    @tail.setter
    def tail(self, tail: str) -> None:
        assert isinstance(tail, str)
        self._tail = tail

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

    def to_json(self):
        data = dict(vars(self))
        for var in ["_tail", "_core", "_leader"]:
            data[var.replace("_", "")] = data[var]
            del data[var]
        data["location"] = str(self.location)
        return data


class Cluster(Feature):
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
        if not self.parent_record:
            raise ValueError("Cluster not contained in record")
        return self.parent_record.get_cluster_number(self)

    def trim_overlapping(self):
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

    def add_cds(self, cds):
        assert isinstance(cds, CDSFeature)
        assert cds.is_contained_by(self), "cds %s outside cluster %s" % (cds, self)
        self.cds_children[cds] = None

    @property
    def cutoff(self):
        return self._cutoff

    @cutoff.setter
    def cutoff(self, cutoff):
        if cutoff is not None:
            cutoff = int(cutoff)
        self._cutoff = cutoff

    @property
    def extent(self):
        return self._extent

    @extent.setter
    def extent(self, extent):
        if extent is not None:
            self._extent = int(extent)

    def get_product_string(self):
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


class NRPSPKSQualifier(list):
    class Domain:
        __slots__ = ["name", "label", "start", "end", "evalue", "bitscore",
                     "predictions"]

        def __init__(self, name, label, start, end, evalue, bitscore):
            self.label = label
            self.name = name
            self.start = start
            self.end = end
            self.evalue = evalue
            self.bitscore = bitscore
            self.predictions = {}

    def __init__(self):
        super().__init__()
        self.type = "uninitialised"
        self.subtypes = []
        self.domains = []
        self.domain_names = []
        self.predictions = {}
        self.cal = 0
        self.at = 0
        self.kr = 0
        self.a = 0
        self.other = 0

    def append(self):
        raise NotImplementedError("Appending to this list won't work, use add_subtype() or add_domain()")

    def extend(self):
        raise NotImplementedError("Extending this list won't work")

    def __len__(self):
        return len(self.subtypes) + len(self.domains)

    def __iter__(self):
        for domain in self.domains:
            yield "NRPS/PKS Domain: %s (%s-%s). E-value: %s. Score: %s;" % (domain.hit_id,
                    domain.query_start, domain.query_end, domain.evalue, domain.bitscore)
        for subtype in self.subtypes:
            yield "NRPS/PKS subtype: %s" % subtype

    def add_subtype(self, subtype):
        assert isinstance(subtype, str)
        self.subtypes.append(subtype)

    def add_domain(self, domain):
        assert not isinstance(domain, str)
        if domain.hit_id == "PKS_AT":
            self.at += 1
            suffix = "_AT%d" % self.at
        elif domain.hit_id == "PKS_KR":
            self.kr += 1
            suffix = "_KR%d" % self.kr
        elif domain.hit_id == "CAL_domain":
            self.cal += 1
            suffix = "_CAL%d" % self.cal
        elif domain.hit_id == "AMP-binding":
            self.a += 1
            suffix = "_A%d" % self.a
        else:
            self.other += 1
            suffix = "_OTHER%d" % self.a

        self.domains.append(NRPSPKSQualifier.Domain(domain.hit_id, suffix,
                domain.query_start, domain.query_end, domain.evalue, domain.bitscore))
        self.domain_names.append(domain.hit_id)


class SecMetQualifier(list):
    def __init__(self, clustertype, domains):
        self.domains = domains  # SecMetResult instance or str
        if domains and not isinstance(domains[0], str):  # SecMetResult
            self.domain_ids = [domain.query_id for domain in self.domains]
        else:  # str
            self.domain_ids = [domain.split()[0] for domain in self.domains]
        self.clustertype = clustertype
        self.kind = "biosynthetic"
        super().__init__()

    def __iter__(self):
        yield "Type: %s" % self.clustertype
        yield "; ".join(map(str, self.domains))
        yield "Kind: %s" % self.kind

    def append(self):
        raise NotImplementedError("Appending to this list won't work")

    def extend(self):
        raise NotImplementedError("Extending this list won't work")

    @staticmethod
    def from_biopython(qualifier):
        domains = None
        clustertype = None
        kind = None
        if len(qualifier) != 3:
            raise ValueError("Cannot parse qualifier: %s" % qualifier)
        for value in qualifier:
            if value.startswith("Type: "):
                clustertype = value.split("Type: ", 1)[0]
            elif value.startswith("Kind: "):
                kind = value.split("Kind: ", 1)[1]
                assert kind == "biosynthetic", kind  # since it's the only kind we have
            else:
                domains = value.split("; ")
        if not domains and clustertype and kind:
            raise ValueError("Cannot parse qualifier: %s" % qualifier)
        return SecMetQualifier(clustertype, domains)

    def __len__(self):
        return 3
