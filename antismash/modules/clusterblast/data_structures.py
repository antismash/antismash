# License: GNU Affero General Public License v3 or later
# A copy of GNU AGPL v3 should have been included in this software package in LICENSE.txt.

from collections import OrderedDict

class ReferenceCluster:
    __slots__ = ["accession", "cluster_label", "proteins", "description",
                 "cluster_type", "tags"]
    def __init__(self, accession, cluster_label, proteins, description,
                 cluster_type, tags):
        self.accession = accession
        self.cluster_label = cluster_label
        self.proteins = proteins
        self.description = description
        self.cluster_type = cluster_type
        self.tags = tags

    def get_name(self):
        return "%s_%s" % (self.accession, self.cluster_label)

    def __getitem__(self, index):
        return [self.proteins, self.description, self.cluster_type, self.tags][index]

    def __str__(self):
        raise ValueError("nope")


class Protein:
    """ Holds details of a protein. Members cannot be added dynamically.
    """
    # At time of writing this class will be instantiated ~7 million times per
    # clusterblast invocation.
    # With those numbers, the memory use of the class without slots: 2.3 Gb
    #                                                and with slots: 0.4 Gb
    __slots__ = ("name", "locus_tag", "location", "strand", "annotations")
    def __init__(self, name, locus_tag, location, strand, annotations):
        self.name = name
        self.locus_tag = locus_tag
        self.location = location
        self.strand = strand
        self.annotations = annotations

    def get_id(self):
        if self.locus_tag and self.locus_tag != "no_locus_tag":
            return self.locus_tag
        return self.name

    def __str__(self):
        if len(self.location.split("-")) != 2:
            raise ValueError("Invalid location in Protein: %s"%self.location)
        tag = self.locus_tag
        if tag == "no_locus_tag":
            tag = self.name
        locations = self.location.replace("-", "\t")
        return "{}\t{}\t{}\t{}\t{}\n".format(tag, self.name, locations,
                                                 self.strand, self.annotations)


class Subject:
    def __init__(self, name, genecluster, start, end, strand, annotation,
                 perc_ident, blastscore, perc_coverage, evalue, locus_tag):
        self.name = name
        self.genecluster = genecluster
        self.start = start
        self.end = end
        self.strand = strand
        self.annotation = annotation
        self.perc_ident = int(perc_ident)
        self.blastscore = int(blastscore)
        self.perc_coverage = float(perc_coverage)
        self.evalue = float(evalue)
        self.locus_tag = locus_tag

    def get_table_string(self):
        return "\t".join([str(i) for i in [self.name, self.perc_ident,
                                           self.blastscore, self.perc_coverage,
                                           self.evalue]])

    @staticmethod
    def from_dict(data):
        args = []
        for key in ["name", "genecluster", "start", "end", "strand",
                    "annotation", "perc_ident", "blastscore", "perc_coverage",
                    "evalue", "locus_tag"]:
            args.append(data[key])
# stop pylint throwing errors because it can't work out the splat operator
# pylint: disable=no-value-for-parameter
        return Subject(*args)
# pylint: enable=no-value-for-parameter


class Query:
    def __init__(self, entry, index):
        parts = entry.split("|")
        self.cluster_number = int(parts[1][1:]) # c1 -> 1
        self.id = parts[4]
        self.entry = entry
        self.subjects = OrderedDict()
        self.cluster_name_to_subjects = {}
        self.index = index

    def add_subject(self, subject):
        self.subjects[subject.name] = subject
        if subject.genecluster not in self.cluster_name_to_subjects:
            self.cluster_name_to_subjects[subject.genecluster] = []
        self.cluster_name_to_subjects[subject.genecluster].append(subject)

    def get_subjects_by_cluster(self, cluster_name):
        return self.cluster_name_to_subjects.get(cluster_name, [])

class Score:
    __slots__ = ("hits", "core_gene_hits", "blast_score", "synteny_score",
                 "core_bonus", "scored_pairings")
    def __init__(self):
        self.hits = 0
        self.core_gene_hits = 0
        self.blast_score = 0
        self.synteny_score = 0
        self.core_bonus = 0
        self.scored_pairings = []

    def sort_score(self):
        """ the pre-existing, unexplained sort weighting """
        if self.core_gene_hits:
            self.core_bonus = 3
        return (self.hits
                + self.core_bonus
                + self.core_gene_hits / 100.
                + self.blast_score / 10e7
                + self.synteny_score)

class MibigEntry:
    def __init__(self, gene_id, gene_description, mibig_cluster,
                mibig_product, percent_id, blast_score, coverage, evalue):
        self.gene_id = gene_id
        self.gene_description = gene_description
        self.mibig_id = mibig_cluster.split("_c")[0]
        self.mibig_product = mibig_product
        self.percent_id = float(percent_id)
        self.blast_score = float(blast_score)
        self.coverage = float(coverage)
        self.evalue = float(evalue)

    @property
    def values(self):
        return [self.gene_id, self.gene_description, self.mibig_id,
                self.mibig_product, self.percent_id, self.blast_score,
                self.coverage, self.evalue]

    def __str__(self):
        return "%s\n" % "\t".join(str(val) for val in self.values)

