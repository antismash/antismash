# License: GNU Affero General Public License v3 or later
# A copy of GNU AGPL v3 should have been included in this software package in LICENSE.txt.

""" A collection of data structures shared by the clusterblast variants """

from collections import OrderedDict
from dataclasses import dataclass
from typing import Any, Dict, List, Self, Tuple, Union

from antismash.common.json import JSONBase, JSONCompatible


@dataclass(slots=True)
class ReferenceCluster:
    """ A reference cluster container, as read from a database of
        antismash-predicted clusters.
    """
    accession: str
    cluster_label: str
    proteins: list[str]
    description: str
    cluster_type: str
    tags: list[str]
    start: int = 0
    end: int = 0

    def __post_init__(self) -> None:
        assert self.cluster_label.startswith("c")
        if "-" not in self.cluster_label:  # no coordinates provided (e.g. subclusterblast)
            start = 0
            end = 0
        else:
            start, end = list(map(int, self.cluster_label.lstrip("c").split("-")))
        self.start = start
        self.end = end

    def __hash__(self) -> int:
        return id(self)

    def get_name(self) -> str:
        """ Returns the name of the cluster, including cluster number """
        return f"{self.accession}_{self.cluster_label}"


@dataclass(slots=True, repr=True)
class Protein(JSONBase):
    """ Holds details of a protein. """
    name: str
    locus_tag: str
    location: str
    strand: str
    annotations: str
    draw_start: int = 0
    draw_end: int = 0

    def __post_init__(self) -> None:
        # set more sane defaults
        try:
            self.draw_start = self.draw_start or self.start
            self.draw_end = self.draw_end or self.end
        except ValueError as err:
            raise ValueError(f"Invalid location in Protein: {self.location}") from err
        assert not self.get_id().startswith("linear")
        assert self.draw_start < self.draw_end

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

    @property
    def start(self) -> int:
        """ The start coordinate of the reference gene within the reference area """
        return int(self.location.split("-")[0])

    @property
    def end(self) -> int:
        """ The end coordinate of the reference gene within the reference area """
        return int(self.location.split("-")[1])

    @classmethod
    def from_json(cls, data: dict[str, Any]) -> Self:
        """ Creates an instance from a raw JSON representation

            Arguments:
                data: the JSON representation to convert

            Returns:
                the newly created instance
        """
        return cls(data["name"], data["locus_tag"], data["location"], data["strand"], data["annotations"],
                   data.get("draw_start", 0), data.get("draw_end", 0))

    def to_json(self) -> JSONCompatible:
        """ Returns a JSON-compatible object with the instance's data """
        return self


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
                 "core_bonus", "scored_pairings", "_similarity",
                 )

    def __init__(self) -> None:
        self.hits = 0
        self.core_gene_hits = 0
        self.blast_score = 0.
        self.synteny_score = 0
        self.core_bonus = 0
        self.scored_pairings: List[Tuple[Query, Subject]] = []
        self._similarity: int = -1

    @property
    def similarity(self) -> int:
        """ The similarity of the pairing, as a percentage in the range 0 to 100 """
        if self._similarity < 0:
            raise ValueError("Similarity not yet set")
        return self._similarity

    @similarity.setter
    def similarity(self, similarity: int) -> None:
        if not 0 <= similarity <= 100:
            raise ValueError(f"Similarity ({similarity}) out of bounds (0-100)")
        self._similarity = similarity

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
