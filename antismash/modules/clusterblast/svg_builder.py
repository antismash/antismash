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
import logging
import os
from typing import Dict, Iterable, List, Set, Tuple, TypeVar

from pysvg.structure import Svg, G as Group
from pysvg.text import Text
from pysvg.builders import ShapeBuilder

from antismash.common.secmet import Feature
from antismash.common.path import get_full_path
from antismash.config import get_config

from .data_structures import Protein, Query, Subject

T = TypeVar('T')


def get_antismash_db_accessions() -> Set[str]:
    """ Find all accession numbers available in the antiSMASH database

        Caches return value for reuse

        Arguments:
            None

        Returns:
            a set of all accessions
    """
    # have we generated them previously
    if not hasattr(get_antismash_db_accessions, "result"):
        filename = get_full_path(__file__, os.path.join('data', 'accessions_in_db.txt'))
        with open(filename, 'r') as handle:
            text = handle.read()

        accessions = text.split('\n')
        if accessions[-1] == '':
            accessions.pop()

        # cache the results
        get_antismash_db_accessions.result = set(accessions)

    # return the cached result
    return get_antismash_db_accessions.result


def generate_distinct_colours(count) -> List[str]:
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
        colours.append("#%s%s%s" % (red, green, blue))
    assert len(colours) == count
    return colours


def sort_groups(query_ids, groups: Iterable[Iterable]) -> List[Iterable]:
    """ Sorts groups into the same order as query_ids. If multiple query_ids are
        in the same group, the earlier id is used for ordering.

        Args:
            query_ids: An iterable of group members defining the sort order.
            groups: The groups to sort, must have a __contains__ method.

        Returns:
            A new list containing the groups in sorted order.
    """

    ordered_groups = []
    found_groups = set()
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


def arrange_colour_groups(query_genes, groups: Iterable[Iterable]) -> List[Iterable]:
    """ Arrange provided groups to be maximally distant from each other.

        Arguments:
            query_genes: the CDS features in the query
            groups: the groupings to rearrange

        Returns:
            a list ordering the original members of groups
    """
    # first sort them
    ordered_groups = sort_groups([gene.get_accession() for gene in query_genes], groups)
    return make_neighbours_distinct(ordered_groups)


def build_colour_groups(query_genes, ranking) -> Dict[str, str]:
    """ Generate a colour for each distinct group of genes.

        A group of genes is distinct if there are no links between any of the
        queries or hits in one group to queries or hits in another group.

        Arguments:
            query_genes: the CDS features from a cluster
            ranking: the cluster-score pairings from a ClusterResult instance

        Returns:
            a dictionary mapping each gene accession to the distinct colour for
            the group the gene belongs to
    """
    # start with a set per query gene with only itself
    groups = {gene.get_accession(): set() for gene in query_genes}
    # populate the sets with the id of any hits matching a query gene
    for _, score in ranking:
        for query, subject in score.scored_pairings:
            groups[query.id].add(query.id)
            groups[query.id].add(subject.name)
    # merge intersecting sets so only disjoint sets remain
    keys = list(groups.keys())
    for i, accession in enumerate(keys):
        groups[accession].add(accession)  # add the query name to it's group
        for other in keys[i + 1:]:
            if groups[accession].intersection(groups[other]):
                groups[other].update(groups[accession])  # merge into later hit
                groups[accession] = groups[other]  # point to merged version

    # only track unique groups with more than one hit, since no match, no colour
    disjoints = set(tuple(group) for group in groups.values() if len(group) > 1)
    # build a reverse lookup from id to group_number
    lookup = {}
    for group, colour in zip(arrange_colour_groups(query_genes, disjoints),
                             generate_distinct_colours(len(disjoints))):
        for name in sorted(group):  # sort to ensure cross-run consistency
            lookup[name] = colour
    return lookup


class Gene:
    """ For constructing the SVG components requried to display a single gene
    """
    def __init__(self, start, end, strand, name, protein=None, product=None):
        self.start = start
        self.end = end
        self.strand = 0
        if not isinstance(strand, int):
            if strand == "+":
                strand = 1
            elif strand == "-":
                strand = -1
        self.strand = strand
        self.name = name
        self.protein = protein
        self.product = product
        self.label = self._get_label().replace("_", " ")
        self.pairings = []
        self.reversed = False

    def _get_label(self):
        if self.protein:
            return self.protein.annotations
        if self.product:
            return self.product
        return self.name

    @staticmethod
    def from_feature(feature) -> 'Gene':  # string because forward reference
        """ Constructs a Gene instance from a CDS feature """
        start = int(feature.location.start)
        end = int(feature.location.end)
        strand = feature.location.strand
        name = feature.get_accession()
        return Gene(start, end, strand, name, product=feature.product)

    @staticmethod
    def from_protein(protein):
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

    def get_end(self):
        """ Returns the maximum edge of the Gene """
        return max([self.start, self.end])

    def _get_description(self):
        """ Returns a HTML fragment describing the Gene """
        description = ['%s[br]Location: %s - %s' % (self.label, self.start, self.end)]
        format_string = ("<br><br><b>BlastP hit with %s</b><br>Percentage identity: %s<br>"
                         "Percentage coverage: %.0f<br>BLAST bit score: %s<br>E-value: %s")
        for query, subject in self.pairings:
            description.append(format_string % (query.id, subject.perc_ident,
                                                float(subject.perc_coverage),
                                                subject.blastscore, subject.evalue))
        return "".join(description)

    def get_arrow_polygon(self, scaling=1., offset=0, base=35, height=10, colour="white"):
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
    def __init__(self, query_cluster_number, ref_cluster_number, accession,
                 description, genes, rank, hits=0, strand=1):
        self.query_cluster_number = query_cluster_number
        self.ref_cluster_number = ref_cluster_number
        self.accession = accession
        self.description = description.replace("_", " ")
        genes = list(genes)
        if isinstance(genes[0], Feature):
            genes = [Gene.from_feature(gene) for gene in genes]
        elif isinstance(genes[0], Protein):
            genes = [Gene.from_protein(gene) for gene in genes]
        else:
            raise TypeError("No conversion from type %s to Gene" % type(genes[0]))
        self.unique_hit_count = hits
        self.genes = sorted(genes, key=lambda gene: gene.start)
        self.overall_strand = strand
        self.reversed = False
        self.start = self.genes[0].get_start()
        self.end = self.genes[-1].get_end()
        self.rank = rank

    @property
    def similarity(self) -> str:
        """ Return a string representing percentage similarity to the query
            cluster
        """
        perc_sim = self.unique_hit_count * 100 / len(self.genes)
        return "%d%% of genes show similarity" % perc_sim

    @property
    def full_description(self) -> str:
        """ Returns a string formatted with accession, cluster number and
            similarity score """
        desc = "%s_%s: %s" % (self.accession, self.ref_cluster_number, self.description)
        if len(desc) > 80:
            desc = desc[:77] + "..."
        return "%s (%s)" % (desc, self.similarity)

    def reverse_strand(self) -> None:
        """ Reverses the entire cluster's directionality, useful when the
            query strand is the opposite direction """
        self.reversed = not self.reversed
        self.overall_strand *= -1
        for gene in self.genes:
            gene.reverse()
        self.genes.reverse()  # to ensure always drawing from left to right

    def __len__(self):
        return abs(self.start - self.end)

    def _add_label(self, group, v_offset):
        acc = Text(self.full_description, 5, 20 + v_offset)
        if self.accession.startswith('BGC'):
            acc = Text('<a xlink:href="http://mibig.secondarymetabolites.org/repository/'
                       + self.accession + '/index.html#cluster-1" target="_blank">'
                       + self.accession + '</a>: '
                       + "%80s (%s)" % (self.description, self.similarity), 5, 20 + v_offset)
        elif self.accession.split("_")[0] in get_antismash_db_accessions():
            acc = Text('<a xlink:href="http://antismash-db.secondarymetabolites.org/output/'
                       + self.accession + '/index.html#cluster-%s" target="_blank">' % self.ref_cluster_number[1:]
                       + self.full_description.replace(":", "</a>:"), 5, 20 + v_offset)
        acc.set_class("clusterblast-acc")
        group.addElement(acc)
        group.setAttribute('label', self.accession)
        group.setAttribute('description', self.full_description)
        group.set_class('clusterblast-cluster')
        return group

    def get_svg_groups(self, h_offset=0, v_offset=0, scaling=1., screenwidth=1024,
                       colours=None, overview=False, prefix="dummy"):
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
            arrow.set_id("%s-%d_%s%d_%s_%s" % (prefix, self.query_cluster_number, label,
                                               self.query_cluster_number, self.rank, i))
            group.addElement(arrow)
            # Can be used for domains
            group.set_id("a%s_00%s" % (self.query_cluster_number, i))
            groups.append(group)
        return groups

    @staticmethod
    def from_reference_cluster(cluster, query_cluster_number, score,
                               reference_proteins, rank, num_hits, strand):
        """ Constructs a Cluster instance from a ReferenceCluster instance """
        proteins = [reference_proteins[protein] for protein in cluster.proteins]
        cluster = Cluster(query_cluster_number, cluster.cluster_label,
                          cluster.accession, cluster.description, proteins,
                          rank, num_hits, strand)
        for query, subject in score.scored_pairings:
            for gene in cluster.genes:
                if gene.name in [subject.name, query.id]:
                    gene.pairings.append((query, subject))
        return cluster


class QueryCluster(Cluster):
    """ A special case of Cluster for the query, with slightly different info
        in the SVG components
    """
    def __init__(self, cluster_feature):
        cluster_number = cluster_feature.get_cluster_number()
        super().__init__(cluster_number, cluster_number,
                         "%s_%d" % (cluster_feature.parent_record.id, cluster_number),
                         "Query sequence",
                         cluster_feature.cds_children, rank=0)

    def get_svg_groups(self, h_offset=0, v_offset=0, scaling=1., screenwidth=1024,
                       colours=None, overview=False, prefix="dummy"):
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
            arrow.set_id("%s-%s_q%s_%s_%s" % (prefix, self.query_cluster_number, index, self.rank, "all"))
            group.addElement(arrow)
        return groups


def determine_strand_of_cluster(cluster, pairings: List[Tuple[Query, Subject]]) -> int:
    """ Determines the strand of a cluster relative to the query cluster.
        Calculated by using the median strand of all linked genes.

        In the case of tie in the counts (e.g. 4 of -1 and 4 of 1) the strand
        of the largest linked gene is used.

        There may also be a tie for the size of largest, so the first will be
        used.

        Arguments:
            cluster: the secmet.Cluster feature to operate on
            pairings: a list of Query,Subject pairs for determining matches

        Returns:
            an int in -1, 0, 1
    """
    name_to_feature = {}
    for cds in cluster.cds_children:
        name_to_feature[cds.get_accession()] = cds

    counted = set()
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
    def __init__(self, cluster_feature, ranking, reference_proteins, prefix):
        if ranking:
            assert reference_proteins
        self.prefix = prefix
        self.query_cluster = QueryCluster(cluster_feature)
        query_cluster_number = cluster_feature.get_cluster_number()
        cluster_limit = get_config().cb_nclusters
        self.colour_lookup = build_colour_groups(cluster_feature.cds_children, ranking[:cluster_limit])
        self.hits = []
        record_prefix = cluster_feature.parent_record.id.split(".", 1)[0]
        num_added = 0
        queries = set()

        for cluster, score in ranking:
            if record_prefix == cluster.accession.split("_", 1)[0]:
                continue
            # determine overall strand direction of hits
            hit_genes = set()
            strand = determine_strand_of_cluster(cluster_feature, score.scored_pairings)
            for query, subject in score.scored_pairings:
                queries.add(query.id)
                hit_genes.add(subject.name)
            cluster = Cluster.from_reference_cluster(cluster, query_cluster_number,
                                                     score, reference_proteins,
                                                     num_added + 1, len(hit_genes),
                                                     strand)
            self.hits.append(cluster)
            num_added += 1
            # obey the cluster display limit from options
            if num_added >= cluster_limit:
                break

        self.max_length = self._size_of_largest_cluster()
        self._organise_strands()

    def get_cluster_descriptions(self) -> List[str]:
        """ Returns all hit cluster descriptions """
        return [cluster.full_description for cluster in self.hits]

    def get_cluster_accessions(self) -> List[str]:
        """ Returns all hit cluster accessions """
        return [cluster.accession for cluster in self.hits]

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

    def get_overview_contents(self, width, height):
        """ Generate an SVG comparing the query to all hit cluster

            Arguments:
                width: the width of the SVG
                height: the height of the SVG

            Returns:
                a string containing the SVG XML
        """
        svg = Svg(x=0, y=0, width=width, height=height)
        viewbox = "0 0 %d %d" % (width, height)
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

    def get_pairing_contents(self, index, width, height) -> str:
        """ Generate an SVG comparing the query to a single hit cluster

            Arguments:
                index: the index of the hit
                width: the width of the SVG
                height: the height of the SVG

            Returns:
                a string containing the SVG XML
        """
        svg = Svg(x=0, y=0, width=width, height=height)
        viewbox = "0 0 %d %d" % (width, height)
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
