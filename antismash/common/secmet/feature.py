# License: GNU Affero General Public License v3 or later
# A copy of GNU AGPL v3 should have been included in this software package in LICENSE.txt.

from collections import OrderedDict
import logging

from Bio.SeqFeature import SeqFeature, FeatureLocation, CompoundLocation

class Feature:
    __slots__ = ["location", "notes", "type", "_qualifiers"]
    def __init__(self, location, feature_type):
        assert isinstance(location, (FeatureLocation, CompoundLocation))
        self.location = location
        self.notes = []
        assert feature_type
        self.type = str(feature_type)
        self._qualifiers = {}

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

    def overlaps_with(self, other):
        assert isinstance(other, Feature)
        return self.location.start in other.location \
               or self.location.end in other.location \
               or other.location.start in self.location \
               or other.location.end in self.location

    def to_biopython(self, qualifiers=None):
        feature = SeqFeature(self.location, type=self.type)
        feature.qualifiers.update(self._qualifiers)
        notes = self._qualifiers.get("note", []) + self.notes
        if qualifiers:
            notes += qualifiers.pop("note", [])
            feature.qualifiers.update(qualifiers)
        if notes:
             # sorting helps with consistency and comparison testing
            feature.qualifiers["note"] = sorted(notes)
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
        if bio_feature.type == "STS" and 2090040 in bio_feature.location:
            assert "db_xref" in bio_feature.qualifiers, bio_feature.location
        if feature is None:
            feature = Feature(bio_feature.location, bio_feature.type)
            if not leftovers:
                leftovers = dict(bio_feature.qualifiers)
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
        self._evalue = None #float
        self.label = None
        self.locus_tag = None
        self._score = None #float

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
    def value(self, evalue):
        self._evalue = float(evalue)

    def to_biopython(self, qualifiers=None):
        mine = OrderedDict()
        if self.label:
            mine["label"] = [self.label]
        if self.score is not None:
            mine["score"] = [str(self.score)]
        if self.evalue is not None:
            mine["evalue"] = [str(self.evalue)]
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
        self.domain_id = None

    def to_biopython(self, qualifiers=None):
        mine = OrderedDict()
        if self.tool:
            mine["aSTool"] = [self.tool]
        if self.domain:
            mine["aSDomain"] = [self.domain]
        if self.domain_id:
            mine["domain_id"] = [self.domain_id]
        if qualifiers:
            mine.update(qualifiers)
        return super().to_biopython(mine)

class CDSMotif(Domain):
    __slots__ = ["motif"]
    def __init__(self, location):
        super().__init__(location, feature_type="CDS_motif")
        self.motif = None # keep

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
    def __init__(self, location, description):
        super().__init__(location, feature_type="PFAM_domain")
        assert isinstance(description, list)
        self.description = description
        self.db_xref = None

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
        feature.db_xref = leftovers.pop("db_xref", None)

        # grab parent optional qualifiers
        super(PFAMDomain, feature).from_biopython(bio_feature, feature=feature, leftovers=leftovers)

        return feature


class AntismashDomain(Domain):
    __slots__ = ["domain_subtype", "specificity"]
    def __init__(self, location):
        super().__init__(location, feature_type="aSDomain")
        self.domain_subtype = None #keep
        self.specificity = None #keep

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

class CDSFeature(Feature):
    __slots__ = ["_translation", "protein_id", "locus_tag", "gene", "product",
                 "transl_table", "sec_met", "aSProdPred", "cluster"]
    _counter = 0
    def __init__(self, location, translation, locus_tag=None, protein_id=None,
                 product=None, gene=None):
        super().__init__(location, feature_type="CDS")
        # mandatory
        #  codon_start
        #  db_xref
        self._translation = str(translation)

        # semi-optional
        self.protein_id = protein_id
        self.locus_tag = locus_tag
        self.gene = gene

        # optional
        self.product = product
        if self.product:
            assert product[0] == "N"
        self.transl_table = None
        self.sec_met = None #SecMetQualifier()
        self.aSProdPred = [] # TODO: shift into nrps sub section?

        if not (protein_id or locus_tag or gene):
            raise ValueError("CDSFeature requires at least one of: gene, protein_id, locus_tag")

        # TODO: add when active site finder completed
        #self.aSASF_choice = []
        #self.aSASF_note = []
        #self.aSASF_prediction = []
        #self.aSASF_scaffold = []

        # runtime-only data
        self.cluster = None

    @property
    def translation(self):
        return self._translation

    @translation.setter
    def translation(self, translation):
        self._translation = str(translation)

    def get_accession(self):
        for val in [self.protein_id, self.gene, self.locus_tag]:
            if val:
                return val
        raise ValueError("%s altered to contain no identifiers" % self)

    def get_name(self):
        for val in [self.locus_tag, self.gene, self.protein_id]:
            if val:
                return val
        raise ValueError("%s altered to contain no identifiers" % self)

    @staticmethod
    def from_biopython(bio_feature, feature=None, leftovers=None):
        if leftovers is None:
            leftovers = Feature.make_qualifiers_copy(bio_feature)
        # grab mandatory qualifiers and create the class
        translation = leftovers.pop("translation")[0]

        # semi-optional qualifiers
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
        feature.sec_met = leftovers.pop("sec_met", None)

        # grab parent optional qualifiers
        super(CDSFeature, feature).from_biopython(bio_feature, feature=feature, leftovers=leftovers)

        return feature

    def to_biopython(self, qualifiers=None):
        mine = OrderedDict()
        # mandatory
        mine["translation"] = [self.translation]
        # optional
        for attr in ["gene", "aSProdPred", "transl_table", "locus_tag",
                     "protein_id", "product"]:
            val = getattr(self, attr)
            if val:
                mine[attr] = [str(val)]
        # since it's already a list
        if self.sec_met:
            mine["sec_met"] = self.sec_met
        # respect qualifiers given to us
        if qualifiers:
            mine.update(qualifiers)
        return super().to_biopython(mine)


class Cluster(Feature):
    __slots__ = ["_extent", "_cutoff", "products", "contig_edge",
                 "detection_rules", "smiles_structure", "probability",
                 "clusterblast", "knownclusterblast", "subclusterblast",
                 "parent_record", "cds_children", "borders"]
    def __init__(self, location, cutoff, extent, products):
        super().__init__(location, feature_type="cluster")

        self._extent = int(extent)
        self._cutoff = int(cutoff)
        assert isinstance(products, list)
        self.products = products

        self.contig_edge = None # hmm_detection borderpredict
        self.detection_rules = []
        self.smiles_structure = None # SMILES string
        self.probability = None  # clusterfinder probability # TODO: unify with notes version

        self.clusterblast = None
        self.knownclusterblast = None
        self.subclusterblast = None

        # for runtime management
        self.parent_record = None
        self.cds_children = []
        self.borders = []

    def get_cluster_number(self):
        if not self.parent_record:
            return None
        return self.parent_record.get_cluster_number(self)

    def add_cds(self, cds):
        assert isinstance(cds, CDSFeature)
        assert cds.location.start in self.location \
            or cds.location.end in self.location
        self.cds_children.append(cds)

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
            text = text.split(":", 1)[1] # strip the leadin
            text = text.split(";") # separate rules
            for rule in text:
                if not rule:
                    continue
                assert ": " in rule, rule
                product, rule = rule.split(": ")
                rule = rule[1:-1] # strip ( )
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
