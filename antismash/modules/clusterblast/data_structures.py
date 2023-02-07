# License: GNU Affero General Public License v3 or later
# A copy of GNU AGPL v3 should have been included in this software package in LICENSE.txt.

""" A collection of data structures shared by the clusterblast variants """

from collections import OrderedDict
from typing import Dict, List, Tuple, Union


class ReferenceCluster:
    """ A reference cluster container, as read from a database of
        antismash-predicted clusters.
    """
    __slots__ = ["accession", "cluster_label", "proteins", "description",
                 "cluster_type", "tags"]

    def __init__(self, accession: str, cluster_label: str, proteins: List[str],
                 description: str, cluster_type: str, tags: List[str]) -> None:
        self.accession = accession
        assert cluster_label.startswith('c')
        self.cluster_label = cluster_label
        self.proteins = proteins
        self.description = description.replace('biosynthetic_gene_cluster', '')
        self.cluster_type = cluster_type
        self.tags = tags

    def get_name(self) -> str:
        """ Returns the name of the cluster, including cluster number """
        return f"{self.accession}_{self.cluster_label}"


class Protein:
    """ Holds details of a protein. Members cannot be added dynamically.
    """
    # At time of writing this class will be instantiated ~7 million times per
    # clusterblast invocation.
    # With those numbers, the memory use of the class without slots: 2.3 Gb
    #                                                and with slots: 0.4 Gb
    __slots__ = ("name", "locus_tag", "location", "strand", "annotations")

    def __init__(self, name: str, locus_tag: str, location: str, strand: str, annotations: str) -> None:
        self.name = name
        self.locus_tag = locus_tag
        self.location = location
        self.strand = strand
        self.annotations = annotations

    def get_id(self) -> str:
        """ Returns best identifier for a Protein """
        if self.locus_tag and self.locus_tag != "no_locus_tag":
            return self.locus_tag
        return self.name

    def __str__(self) -> str:
        if len(self.location.split("-")) != 2:
            raise ValueError(f"Invalid location in Protein: {self.location}")
        tag = self.locus_tag
        if tag == "no_locus_tag":
            tag = self.name
        locations = self.location.replace("-", "\t")
        return f"{tag}\t{self.name}\t{locations}\t{self.strand}\t{self.annotations}\n"


class Subject:
    """ Holds details of a subject as reported by BLAST """
    def __init__(self, name: str, genecluster: str, start: int, end: int, strand: str, annotation: str,
                 perc_ident: int, blastscore: int, perc_coverage: float, evalue: float, locus_tag: str) -> None:
        self.name = name
        self.genecluster = genecluster
        self.start = int(start)
        self.end = int(end)
        self.strand = str(strand)
        self.annotation = annotation
        self.perc_ident = int(perc_ident)
        self.blastscore = int(blastscore)
        self.perc_coverage = float(perc_coverage)
        self.evalue = float(evalue)
        self.locus_tag = str(locus_tag)

    def __len__(self) -> int:
        return abs(int(self.start) - int(self.end))

    def get_table_string(self) -> str:
        """ Returns a string of the Subject suitable for use in writing text tables """
        return "\t".join([str(i) for i in [self.name, self.perc_ident,
                                           self.blastscore, self.perc_coverage,
                                           self.evalue]])

    @staticmethod
    def from_dict(data: Dict[str, Union[str, int, float]]) -> "Subject":
        """ Recreates a Subject instance from a JSON formatted subject """
        args = []
        for key in ["name", "genecluster", "start", "end", "strand",
                    "annotation", "perc_ident", "blastscore", "perc_coverage",
                    "evalue", "locus_tag"]:
            args.append(data[key])
        return Subject(*args)  # type: ignore # pylint: disable=no-value-for-parameter


class Query:
    """ Holds details of a query as reported by blast, with links to subjects
        that were connected to the query
    """
    def __init__(self, entry: str, index: int) -> None:
        parts = entry.split("|")
        self.cluster_number = int(parts[1][1:])  # c1 -> 1
        self.id = parts[4]  # accession
        self.entry = entry
        self.subjects: Dict[str, Subject] = OrderedDict()
        self.cluster_name_to_subjects: Dict[str, List[Subject]] = {}
        self.index = index

    def add_subject(self, subject: Subject) -> None:
        """ Adds a Subject to the Query linkings """
        self.subjects[subject.name] = subject
        if subject.genecluster not in self.cluster_name_to_subjects:
            self.cluster_name_to_subjects[subject.genecluster] = []
        self.cluster_name_to_subjects[subject.genecluster].append(subject)

    def get_subjects_by_cluster(self, cluster_name: str) -> List[Subject]:
        """ Returns a list of Subjects that shared the cluster name """
        return self.cluster_name_to_subjects.get(cluster_name, [])


class Score:
    """ A multi-part score for a cluster """
    __slots__ = ("hits", "core_gene_hits", "blast_score", "synteny_score",
                 "core_bonus", "scored_pairings")

    def __init__(self) -> None:
        self.hits = 0
        self.core_gene_hits = 0
        self.blast_score = 0.
        self.synteny_score = 0
        self.core_bonus = 0
        self.scored_pairings: List[Tuple[Query, Subject]] = []

    @property
    def score(self) -> int:
        """ The cluster's score """
        if self.core_gene_hits:
            self.core_bonus = 3
        return self.hits + self.core_bonus + self.core_gene_hits + self.synteny_score

    def sort_score(self) -> Tuple[int, float]:
        """ For sorting purposes, sort first by score and solve any ties by
            cumulative blast score
        """
        return (self.score, self.blast_score)


class MibigEntry:
    """ A container for tracking similarity to a MIBiG entry """
    def __init__(self, gene_id: str, gene_description: str, mibig_id: str, mibig_cluster_number: int,
                 mibig_product: str, percent_id: float, blast_score: float,
                 coverage: float, evalue: float) -> None:
        self.gene_id = gene_id
        self.gene_description = gene_description
        self.mibig_id = mibig_id
        self.mibig_cluster_number = mibig_cluster_number
        self.mibig_product = mibig_product
        self.percent_id = float(percent_id)
        self.blast_score = float(blast_score)
        self.coverage = float(coverage)
        self.evalue = float(evalue)

    @property
    def values(self) -> List[Union[str, float]]:
        """ a list of all class member values in constructor arg order, for
            simplifying conversion to and from JSON
        """
        return [self.gene_id, self.gene_description, self.mibig_id, self.mibig_cluster_number,
                self.mibig_product, self.percent_id, self.blast_score,
                self.coverage, self.evalue]

    def __str__(self) -> str:
        return "%s\n" % "\t".join(str(val) for val in self.values)
