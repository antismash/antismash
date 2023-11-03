# License: GNU Affero General Public License v3 or later
# A copy of GNU AGPL v3 should have been included in this software package in LICENSE.txt.

""" Responsible for the construction of Clusterblast related SVGs.

    Cluster query genes that were matched to any of the genes in a reference
    cluster will be coloured along with the reference genes that matched.
    If a gene would be given two colours, all genes sharing that colour are
    given a single grouped colour, so multiple query genes and reference genes
    may have the same colour.

    Colours are generated to be consistent across multiple runs, however if the
    number of distinct groups changes, the colours will also change.

    Colours for neighbouring genes in the query should be as distant as
    possible, but with enough hits the distance may end up small anyway.

    Information is stored for tooltips, including details of the gene itself
    and, for reference genes, which genes in the query resulted in the hit.


    In order to construct all SVGs, only one instance of ClusterSVGBuilder is
    required for each group of related clusters.
"""

import colorsys
from typing import Dict, Iterable, List, Set, Tuple, TypeVar, Union

from pysvg.builders import ShapeBuilder
from pysvg.shape import Polygon
from pysvg.structure import Svg, G as Group
from pysvg.text import Text

from antismash.common import secmet
from antismash.config import get_config

from .data_structures import Protein, Query, Subject, ReferenceCluster, Score

T = TypeVar('T')


def generate_distinct_colours(count: int) -> List[str]:
    """ Generates `count` non-white colours that are as removed from
        each other as possible while using the same saturation and colour values

        Arguments:
            count: the number of non-white colours required

        Returns:
            a list of colour strings in # hex format
    """
    rgbs = [colorsys.hsv_to_rgb(i/count, .9, .85) for i in range(count)]
    colours = []
    for rgb in rgbs:
        red, green, blue = [str(hex(int(i*255)))[2:] for i in rgb]  # [2:] strips the 0x
        colours.append(f"#{red}{green}{blue}")
    assert len(colours) == count
    return colours


def sort_groups(query_ids: Iterable[T], groups: Set[Tuple[T, ...]]) -> List[Tuple[T, ...]]:
    """ Sorts groups into the same order as query_ids. If multiple query_ids are
        in the same group, the earlier id is used for ordering.

        Args:
            query_ids: An iterable of group members defining the sort order.
            groups: The groups to sort

        Returns:
            A new list containing the groups in sorted order.
    """

    ordered_groups = []
    found_groups: Set[int] = set()
    for query_id in query_ids:
        for group in groups:
            if query_id in group:
                if id(group) not in found_groups:
                    ordered_groups.append(group)
                    found_groups.add(id(group))  # id since sets are unhashable
    return ordered_groups


def make_neighbours_distinct(groups: List[T]) -> List[T]:
    """ Rearranges the incoming list such that no neighbours of the original
        are neighbours in the result. E.g. [0,1,2,3,4] -> [0,2,4,1,3]

        Arguments:
            groups: a collection of values to rearrange to be distant

        Returns:
             a list containing the members of the original container, with each
             member being as distant from it's original neighbours as possible
    """
    spaced_groups = []
    for i in range(4):
        spaced_groups.extend(groups[i::4])
    return spaced_groups


def arrange_colour_groups(accessions: List[str],
                          groups: Set[Tuple[str, ...]]) -> List[Tuple[str, ...]]:
    """ Arrange provided groups to be maximally distant from each other.

        Arguments:
            accessions: the names of features in the query
            groups: the groupings to rearrange

        Returns:
            a list ordering the original members of groups
    """
    # first sort them
    ordered_groups = sort_groups(accessions, groups)
    return make_neighbours_distinct(ordered_groups)


def build_colour_groups(query_cds_features: List[secmet.CDSFeature],
                        ranking: List[Tuple[ReferenceCluster, Score]]) -> Dict[str, str]:
    """ Generate a colour for each distinct group of genes.

        A group of genes is distinct if there are no links between any of the
        queries or hits in one group to queries or hits in another group.

        Arguments:
            query_cds_features: the CDS features from a cluster
            ranking: the cluster-score pairings from a ClusterResult instance

        Returns:
            a dictionary mapping each gene accession to the distinct colour for
            the group the gene belongs to
    """
    # start with a set per query gene with only itself
    accessions = [cds.get_name() + " QUERY" for cds in query_cds_features]
    groups: Dict[str, Set[str]] = {accession: set() for accession in accessions}
    # populate the sets with the id of any hits matching a query gene
    for _, score in ranking:
        for query, subject in score.scored_pairings:
            q_id = query.id + " QUERY"
            groups[q_id].add(query.id)
            groups[q_id].add(subject.name)
    # merge intersecting sets so only disjoint sets remain
    for i, accession in enumerate(accessions):
        groups[accession].add(accession)  # add the query name to it's group
        for other in accessions[i + 1:]:
            if groups[accession].intersection(groups[other]):
                groups[other].update(groups[accession])  # merge into later hit
                groups[accession] = groups[other]  # point to merged version

    # only track unique groups with more than one hit, since no match, no colour
    disjoints = set(tuple(group) for group in groups.values() if len(group) > 1)
    # build a reverse lookup from id to group_number
    lookup = {}
    for group, colour in zip(arrange_colour_groups(accessions, disjoints),
                             generate_distinct_colours(len(disjoints))):
        for name in sorted(group):  # sort to ensure cross-run consistency
            if name.endswith(" QUERY"):
                name = name[:-6]
            lookup[name] = colour
    return lookup


class Gene:
    """ For constructing the SVG components requried to display a single gene
    """
    def __init__(self, start: int, end: int, strand: Union[int, str], name: str,
                 protein: Protein = None, product: str = None) -> None:
        self.start = start
        self.end = end
        self.strand = 0
        if isinstance(strand, str):
            if strand == "+":
                self.strand = 1
            elif strand == "-":
                self.strand = -1
        else:
            assert isinstance(strand, int), type(strand)
            self.strand = strand
        self.name = name
        self.protein = protein
        self.product = product
        self.label = self._get_label().replace("_", " ")
        self.pairings: List[Tuple[Query, Subject]] = []
        self.reversed = False

    def _get_label(self) -> str:
        if self.protein:
            return self.protein.annotations
        if self.product:
            return self.product
        return self.name

    @staticmethod
    def from_feature(feature: secmet.CDSFeature) -> 'Gene':  # string because forward reference
        """ Constructs a Gene instance from a CDS feature """
        start = int(feature.location.start)
        end = int(feature.location.end)
        strand = feature.location.strand
        name = feature.get_name()
        return Gene(start, end, strand, name, product=feature.product)

    @staticmethod
    def from_protein(protein: Protein) -> "Gene":
        """ Constructs a Gene instance from a Protein instance """
        strand = protein.strand
        name = protein.get_id()
        start, end = [int(i) for i in protein.location.split("-")]
        return Gene(start, end, strand, name, protein=protein)

    def reverse(self) -> None:
        """ Flips the strand and start/end locations of the Gene """
        self.reversed = not self.reversed
        self.start, self.end = self.end, self.start
        self.strand *= -1

    def get_start(self) -> int:
        """ Returns the minimum edge of the Gene """
        return min([self.start, self.end])

    def get_end(self) -> int:
        """ Returns the maximum edge of the Gene """
        return max([self.start, self.end])

    def _get_description(self) -> str:
        """ Returns a HTML fragment describing the Gene """
        description = [f"<b>{self.name}</b><br>{self.label}<br>Location: {self.start} - {self.end}"]
        format_string = ("<br><br><b>BlastP hit with %s</b><br>Percentage identity: %s<br>"
                         "Percentage coverage: %.0f<br>BLAST bit score: %s<br>E-value: %s")
        for query, subject in self.pairings:
            description.append(format_string % (query.id, subject.perc_ident,
                                                float(subject.perc_coverage),
                                                subject.blastscore, subject.evalue))
        return "".join(description)

    def get_arrow_polygon(self, scaling: float = 1., offset: int = 0, base: int = 35,
                          height: int = 10, colour: str = "white") -> Polygon:
        """ Returns an SVG polygon shaped like an arrow that represents the Gene """
        if self.reversed:
            start = int((offset - self.get_start()) * scaling)
            end = int((offset - self.get_end()) * scaling)
        else:
            start = int((self.get_start() + offset) * scaling)
            end = int((self.get_end() + offset) * scaling)
        start, end = sorted([start, end])

        builder = ShapeBuilder()
        arrow_size = height // 2
        if abs(start - end) < arrow_size:
            if self.strand > -1:
                points = [(start, base),
                          (end, base - arrow_size),
                          (start, base - height),
                          (start, base)]
            else:
                points = [(start, base - arrow_size),
                          (end, base - height),
                          (end, base),
                          (start, base - arrow_size)]
        else:
            if self.strand > -1:
                arrowstart = end - arrow_size
                points = [(start, base),
                          (arrowstart, base),
                          (end, base - arrow_size),
                          (arrowstart, base - height),
                          (start, base - height),
                          (start, base)]
            else:
                arrowstart = start + arrow_size
                points = [(start, base - arrow_size),
                          (arrowstart, base - height),
                          (end, base - height),
                          (end, base),
                          (arrowstart, base),
                          (start, base - arrow_size)]
        arrow = builder.createPolygon(strokewidth=1, stroke='black', fill=colour,
                                      points=builder.convertTupleArrayToPoints(points))
        locus_tag = self.name
        if self.protein:
            locus_tag = self.protein.get_id()
        arrow.setAttribute('description', self._get_description())
        arrow.setAttribute('locus_tag', locus_tag)
        arrow.set_class('clusterblast-orf')
        return arrow


class Cluster:
    """ For constructing all SVG components required to represent a cluster
    """
    def __init__(self, region_number: int, ref_cluster_number: str, accession: str,
                 description: str, features: Union[List[Protein], List[secmet.CDSFeature]], rank: int,
                 cluster_type: str, hits: int = 0, strand: int = 1, prefix: str = "general") -> None:
        self.region_number = region_number
        self.ref_cluster_number = ref_cluster_number.lstrip('c').replace("<", "").replace(">", "")
        self.accession = accession
        self.description = description.replace("_", " ")
        self.cluster_type = cluster_type
        if isinstance(features[0], secmet.CDSFeature):
            genes = []
            for feature in features:
                assert isinstance(feature, secmet.CDSFeature), type(feature)
                genes.append(Gene.from_feature(feature))
        elif isinstance(features[0], Protein):
            genes = []
            for feature in features:
                assert isinstance(feature, Protein), type(feature)
                genes.append(Gene.from_protein(feature))
        else:
            raise TypeError(f"No conversion from type {type(features[0])} to Gene")
        self.unique_hit_count = int(hits)
        self.genes: List[Gene] = sorted(genes, key=lambda gene: gene.start)
        self.overall_strand = strand
        self.reversed = False
        self.start = int(self.genes[0].get_start())
        self.end = int(self.genes[-1].get_end())
        self.rank = rank
        self.prefix = prefix

    @property
    def similarity(self) -> int:
        """ Returns the percentage similarity to the query cluster as an int 0-100
        """
        return self.unique_hit_count * 100 // len(self.genes)

    @property
    def similarity_string(self) -> str:
        """ Return a string representing percentage similarity to the query
            cluster
        """
        return f"{self.similarity}% of genes show similarity"

    @property
    def full_description(self) -> str:
        """ Returns a string formatted with accession, cluster number and
            similarity score """
        if self.prefix == "general":
            desc = f"{self.accession} ({self.ref_cluster_number}): {self.description}"
        else:
            desc = f"{self.accession}_c{self.ref_cluster_number}: {self.description}"
        if len(desc) > 80:
            desc = desc[:77] + "..."
        return f"{desc} ({self.similarity_string}), {self.cluster_type}"

    def reverse_strand(self) -> None:
        """ Reverses the entire cluster's directionality, useful when the
            query strand is the opposite direction """
        self.reversed = not self.reversed
        self.overall_strand *= -1
        for gene in self.genes:
            gene.reverse()
        self.genes.reverse()  # to ensure always drawing from left to right

    def __len__(self) -> int:
        return abs(self.start - self.end)

    def _add_label(self, group: Group, v_offset: int) -> Group:
        acc = Text(self.full_description, 5, 20 + v_offset)
        if self.prefix == "knownclusterblast":
            desc = f"{self.description:80} ({self.similarity_string}), {self.cluster_type}"
            acc = Text((
                '<a xlink:href="https://mibig.secondarymetabolites.org/go/'
                f'{self.accession}/{self.ref_cluster_number}"'
                ' target="_blank">'
                f"{self.accession}</a>: {desc}"
            ), 5, 20 + v_offset)
        elif self.prefix == "general":
            start, end = self.ref_cluster_number.split("-")
            query = f"record={self.accession}&start={start}&end={end}".replace("&", "&amp;")
            acc = Text(f'<a xlink:href="https://antismash-db.secondarymetabolites.org/area?{query}'
                       + '" target="_blank">'
                       + self.full_description.replace(":", "</a>:", 1), 5, 20 + v_offset)
        # Don't do any linking for subclusterblast

        acc.set_class("clusterblast-acc")
        group.addElement(acc)
        group.setAttribute('label', self.accession)
        group.setAttribute('description', self.full_description)
        group.set_class('clusterblast-cluster')
        return group

    def get_svg_groups(self, h_offset: int = 0, v_offset: int = 0, scaling: float = 1.,
                       screenwidth: int = 1024, colours: Dict[str, str] = None,
                       overview: bool = False, prefix: str = "dummy") -> List[Group]:
        """ Returns all SVG elements required to draw the Cluster """
        if not colours:
            colours = {}
        groups = []
        group = Group()
        self._add_label(group, v_offset)
        line_y = 35 + v_offset
        group.addElement(ShapeBuilder().createLine(10, line_y,
                                                   10 + (screenwidth * 0.75), line_y,
                                                   strokewidth=1, stroke="grey"))
        groups.append(group)
        # Add gene arrows
        arrow_y = line_y + 5
        offset = h_offset - self.start
        if self.reversed:
            offset = h_offset + self.end
        for i, gene in enumerate(self.genes):
            group = Group()
            arrow = gene.get_arrow_polygon(scaling=scaling, offset=offset,
                                           base=arrow_y,
                                           colour=colours.get(gene.name, "white"))
            if overview:
                label = "all_"
            else:
                label = "h"
            arrow.set_id(f"{prefix}-{self.region_number}_{label}{self.region_number}_{self.rank}_{i}")
            group.addElement(arrow)
            # Can be used for domains
            group.set_id(f"a{self.region_number}_00{i}")
            groups.append(group)
        return groups

    @staticmethod
    def from_reference_cluster(cluster: ReferenceCluster, region_number: int, score: Score,
                               reference_proteins: Dict[str, Protein], rank: int, num_hits: int,
                               strand: int, prefix: str) -> "Cluster":
        """ Constructs a Cluster instance from a ReferenceCluster instance """
        proteins = [reference_proteins[protein] for protein in cluster.tags]
        svg_cluster = Cluster(region_number, str(cluster.cluster_label),
                              cluster.accession, cluster.description, proteins,
                              rank, cluster.cluster_type, num_hits, strand, prefix)
        for query, subject in score.scored_pairings:
            for gene in svg_cluster.genes:
                if gene.name in [subject.name, query.id]:
                    gene.pairings.append((query, subject))
        return svg_cluster


class QueryRegion(Cluster):
    """ A special case of Cluster for the query region, with slightly different
        info in the SVG components
    """
    def __init__(self, region_feature: secmet.Region) -> None:
        region_number = region_feature.get_region_number()
        super().__init__(region_number, str(region_number),
                         f"{region_feature.parent_record.id}_{region_number}",
                         "Query sequence",
                         list(region_feature.cds_children), rank=0, cluster_type="query")

    def get_svg_groups(self, h_offset: int = 0, v_offset: int = 0, scaling: float = 1.,
                       screenwidth: int = 1024, colours: Dict[str, str] = None,
                       overview: bool = False, prefix: str = "dummy") -> List[Group]:
        """ Returns all SVG elements required to draw the Cluster """
        if not colours:
            colours = {}
        groups = []
        group = Group()
        acc = Text(self.description, 5, 20)
        acc.set_class("clusterblast-acc")
        group.addElement(acc)
        line_y = 35 + v_offset
        group.addElement(ShapeBuilder().createLine(10, line_y, 10 + (screenwidth * 0.75),
                                                   line_y, strokewidth=1, stroke="grey"))
        group.setAttribute('label', self.description)
        group.set_class('clusterblast-cluster')
        groups.append(group)

        base = line_y + 5
        offset = h_offset + 10 - self.start  # 10 for margin
        for index, gene in enumerate(self.genes):
            arrow = gene.get_arrow_polygon(scaling=scaling, offset=offset,
                                           base=base, colour=colours.get(gene.name, "white"))
            arrow.set_id(f"{prefix}-{self.region_number}_q{index}_{self.rank}_all")
            group.addElement(arrow)
        return groups


def determine_strand_of_cluster(region: secmet.Region, pairings: List[Tuple[Query, Subject]]) -> int:
    """ Determines the strand of a cluster relative to the query cluster.
        Calculated by using the median strand of all linked genes.

        In the case of tie in the counts (e.g. 4 of -1 and 4 of 1) the strand
        of the largest linked gene is used.

        There may also be a tie for the size of largest, so the first will be
        used.

        Arguments:
            cluster: the secmet.Region feature to operate on
            pairings: a list of Query,Subject pairs for determining matches

        Returns:
            an int in -1, 0, 1
    """
    name_to_feature = {}
    for cds in region.cds_children:
        name_to_feature[cds.get_name()] = cds

    counted: Set[str] = set()
    strand = 0
    largest_size = 0
    strand_of_largest = 0
    for query, subject in pairings:
        # only count a gene once for purposes of determining strand
        if subject.name not in counted:
            subject_strand = 1
            if subject.strand == "-":
                subject_strand = -1
            if subject_strand == name_to_feature[query.id].location.strand:
                strand += 1
            else:
                strand -= 1
            if len(subject) > largest_size and subject_strand != 0:
                largest_size = len(subject)
                strand_of_largest = subject_strand
        counted.add(subject.name)
    if strand < 0:
        strand = -1
    elif strand > 0:
        strand = 1
    else:  # a tiebreaker is required
        strand = strand_of_largest
    return strand


class ClusterSVGBuilder:
    """ Constructs SVGs for both query cluster and matching clusters, in both
        pairwise and combined forms
    """
    def __init__(self, region: secmet.Region, ranking: List[Tuple[ReferenceCluster, Score]],
                 reference_proteins: Dict[str, Protein], prefix: str) -> None:
        if ranking:
            assert reference_proteins
        self.prefix = prefix
        self.query_cluster = QueryRegion(region)
        region_number = region.get_region_number()
        cluster_limit = get_config().cb_nclusters
        assert len(ranking) <= cluster_limit, f"{len(ranking)} <= {cluster_limit}"
        self.colour_lookup = build_colour_groups(list(region.cds_children), ranking)
        self.hits: List[Cluster] = []
        num_added = 0
        queries = set()
        for cluster, score in ranking:
            # determine overall strand direction of hits
            hit_genes = set()
            strand = determine_strand_of_cluster(region, score.scored_pairings)
            for query, subject in score.scored_pairings:
                queries.add(query.id)
                hit_genes.add(subject.name)
            svg_cluster = Cluster.from_reference_cluster(cluster, region_number,
                                                         score, reference_proteins,
                                                         num_added + 1, len(hit_genes),
                                                         strand, self.prefix)
            self.hits.append(svg_cluster)
            num_added += 1

        self.max_length = self._size_of_largest_cluster()
        self._organise_strands()

    def get_cluster_descriptions(self) -> List[str]:
        """ Returns all hit cluster descriptions """
        return [cluster.full_description for cluster in self.hits]

    def get_cluster_accessions(self) -> List[str]:
        """ Returns all hit cluster accessions """
        return [cluster.accession for cluster in self.hits]

    def get_cluster_similarities(self) -> List[int]:
        """ Returns all hit cluster similarity percentages """
        return [cluster.similarity for cluster in self.hits]

    def _organise_strands(self) -> None:
        """ Attempt to make all hits have a uniform overall direction """
        for cluster in self.hits:
            if cluster.overall_strand == -1:  # this is relative to the query
                cluster.reverse_strand()

    def _size_of_largest_cluster(self) -> int:
        query_length = len(self.query_cluster)
        length = query_length
        for cluster in self.hits:
            if len(cluster) > length:
                length = len(cluster)
        min_scale = get_config().cb_min_homology_scale
        # if this would shrink the query too much, use the minimum allowed
        if query_length / length < min_scale:
            length = int(query_length / min_scale)
        return length

    def get_overview_contents(self, width: int, height: int) -> str:
        """ Generate an SVG comparing the query to all hit cluster

            Arguments:
                width: the width of the SVG
                height: the height of the SVG

            Returns:
                a string containing the SVG XML
        """
        svg = Svg(x=0, y=0, width=width, height=height)
        viewbox = f"0 0 {width} {height}"
        svg.set_viewBox(viewbox)
        svg.set_preserveAspectRatio("none")
        scaling = (width - 20) / self.max_length  # -20 for margins
        offset = (self.max_length - len(self.query_cluster)) // 2
        for group in self.query_cluster.get_svg_groups(h_offset=offset, scaling=scaling,
                                                       colours=self.colour_lookup, overview=True,
                                                       prefix=self.prefix):
            svg.addElement(group)
        for index, cluster in enumerate(self.hits):
            for group in cluster.get_svg_groups(v_offset=50 * (index + 1),
                                                h_offset=(self.max_length - len(cluster)) // 2,
                                                scaling=scaling, colours=self.colour_lookup, overview=True,
                                                prefix=self.prefix):
                svg.addElement(group)
        return svg.getXML()

    def get_pairing_contents(self, index: int, width: int, height: int) -> str:
        """ Generate an SVG comparing the query to a single hit cluster

            Arguments:
                index: the index of the hit
                width: the width of the SVG
                height: the height of the SVG

            Returns:
                a string containing the SVG XML
        """
        svg = Svg(x=0, y=0, width=width, height=height)
        viewbox = f"0 0 {width} {height}"
        svg.set_viewBox(viewbox)
        svg.set_preserveAspectRatio("none")
        max_length = max([len(self.query_cluster), len(self.hits[index])])
        scaling = (width - 20) / max_length
        for i, cluster in enumerate([self.query_cluster, self.hits[index]]):
            for group in cluster.get_svg_groups(v_offset=50 * i,
                                                h_offset=(max_length - len(cluster)) // 2,
                                                scaling=scaling, colours=self.colour_lookup):
                svg.addElement(group)
        return svg.getXML()
