# License: GNU Affero General Public License v3 or later
# A copy of GNU AGPL v3 should have been included in this software package in LICENSE.txt.

import colorsys
import logging
import os

from pysvg.structure import Svg, G as Group
from pysvg.text import Text
from pysvg.builders import ShapeBuilder

from antismash.common.secmet import Feature
from antismash.common.path import get_full_path
from antismash.config.args import Config

from .data_structures import Protein

def get_antismash_db_accessions():
    ''' Returns a set of all accession numbers available in the antiSMASH database

        Caches return value for reuse
    '''
    # have we generated them previously
    if not hasattr(get_antismash_db_accessions, "result"):
        filename = get_full_path(__file__, os.path.join('data', 'accessions_in_db.txt'))
        with open(filename, 'r') as handle:
            text = handle.read()

        accessions = text.split('\n')
        if accessions[-1] == '':
            accessions.pop()

        get_antismash_db_accessions.result = set(accessions)

    return get_antismash_db_accessions.result


def generate_distinct_colours(count):
    count += 1 # include white
    rgbs = [colorsys.hsv_to_rgb(i/count, .9, .85) for i in range(count)]
    colours = []
    for rgb in rgbs:
        r, g, b = [str(hex(int(i*255)))[2:] for i in rgb] # [2:] strips the 0x
        colours.append("#%s%s%s" % (r, g, b))
    assert len(colours) >= count
    return colours


def sort_groups(query_ids, groups):
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
                    found_groups.add(id(group)) #since sets are unhashable
    return ordered_groups


def make_neighbours_distinct(groups):
    """ Rearranges the incoming iterable such that no neighbours of the original
        are neighbours in the result. E.g. [0,1,2,3,4] -> [0,2,4,1,3]

        Returns a new list containing the members of the original
    """
    # if there's only 2 groups, we can't fix that, so just return it as is
    if len(groups) < 2:
        return list(groups)
    elif len(groups) < 4: # to avoid removing all groups with division below
        spaced_groups = list(groups)[::2] + list(groups)[1::2]

    spaced_groups = []
    step_size = len(groups) // 4
    for i in range(step_size):
        for group in groups[i::step_size]:
            spaced_groups.append(group)
    return spaced_groups


def arrange_colour_groups(query_genes, groups):
    # first sort them
    ordered_groups = sort_groups([gene.get_accession() for gene in query_genes], groups)
    return make_neighbours_distinct(ordered_groups)


def build_colour_groups(query_genes, ranking):
    # start with a set per query gene with only itself
    groups = {gene.get_accession() : set() for gene in query_genes}
    # populate the sets with the id of any hits matching a query gene
    for _, score in ranking:
        for query, subject in score.scored_pairings:
            groups[query.id].add(subject.name)
    # merge intersecting sets so only disjoint sets remain
    keys = list(groups.keys())
    for i, accession in enumerate(keys):
        groups[accession].add(accession) # add the query name to it's group
        for other in keys[i + 1:]:
            if groups[accession].intersection(groups[other]):
                groups[other].update(groups[accession]) # merge into later hit
                groups[accession] = groups[other] # point to merged version

    # generate a colour for each set, as distinct as possible
    disjoints = set(tuple(group) for group in groups.values() if len(group) > 1) # no match, no colour
    colours = generate_distinct_colours(len(disjoints))
    # build a reverse lookup from id to group_number
    lookup = {}
    for group, colour in zip(arrange_colour_groups(query_genes, disjoints), colours):
        for name in group:
            lookup[name] = colour
    return lookup


class Gene:
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
    def from_feature(feature):
        start = int(feature.location.start)
        end = int(feature.location.end)
        strand = feature.location.strand
        name = feature.get_accession()
        return Gene(start, end, strand, name, product=feature.product)

    @staticmethod
    def from_protein(protein):
        strand = protein.strand
        name = protein.get_id()
        start, end = [int(i) for i in protein.location.split("-")]
        return Gene(start, end, strand, name, protein=protein)

    def reverse(self):
        self.reversed = not self.reversed
        self.start, self.end = self.end, self.start
        self.strand *= -1

    def get_start(self):
        if self.reversed:
            return max([self.start, self.end])
        return min([self.start, self.end])

    def get_end(self):
        if self.reversed:
            return min([self.start, self.end])
        return max([self.start, self.end])

    def get_block_polygon(self, start, end, base=35, height=10, colour="white"):
        builder = ShapeBuilder()
        points = []
        for x in [start, end]:
            for y in [base, base + height]:
                points.append((x, y))
        block = builder.createPolygon(strokewidth=1, stroke='black', fill=colour,
                            points=builder.convertTupleArrayToPoints(points))
        locus_tag = self.name
        if self.protein:
            locus_tag = self.protein.get_id()
        block.setAttribute('description', self._get_description())
        block.setAttribute('locus_tag', locus_tag)
        block.set_class('clusterblast-orf')
        return block

    def _get_description(self):
        description = ['%s[br]Location: %s - %s' % (self.label, self.start, self.end)]
        for query, subject in self.pairings:
            description.append('<br><br><b>BlastP hit with %s</b><br>Percentage identity: %s<br>Percentage coverage: %s<br>BLAST bit score: %s<br>E-value: %s' % (query.id, subject.perc_ident, subject.perc_coverage, subject.blastscore, subject.evalue))
        return "".join(description)

    def get_arrow_polygon(self, scaling=1., offset=0, base=35, height=10, colour="white"):
        start = int((self.get_start() + offset) * scaling)
        end = int((self.get_end() + offset) * scaling)
        if self.reversed:
            start = int((offset - self.get_start()) * scaling)
            end = int((offset - self.get_end()) * scaling)

        if not self.strand:
            return self.get_block_polygon(start, end, height)

        builder = ShapeBuilder()
        arrow_size = height // 2
        if abs(start - end) < arrow_size:
            if self.strand > -1:
                points = [(start, base),
                          (end, base - arrow_size),
                          (start, base - height),
                          (start, base)
                         ]
            else:
                points = [(start, base - arrow_size),
                          (end, base - height),
                          (end, base),
                          (start, base - arrow_size)
                         ]
        else:
            if self.strand > -1:
                arrowstart = end - arrow_size
                points = [(start, base),
                          (arrowstart, base),
                          (end, base - arrow_size),
                          (arrowstart, base - height),
                          (start, base - height),
                          (start, base)
                         ]
            else:
                arrowstart = start + arrow_size
                points = [(start, base - arrow_size),
                          (arrowstart, base - height),
                          (end, base - height),
                          (end, base),
                          (arrowstart, base),
                          (start, base - arrow_size)
                         ]
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
    def __init__(self, query_cluster_number, ref_cluster_number, accession, description, genes, rank, hits=0, strand=1):
        self.query_cluster_number = query_cluster_number
        self.ref_cluster_number = ref_cluster_number
        self.accession = accession
        self.description = description.replace("_", " ")
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
    def full_description(self):
        desc = self.description
        if len(desc) > 80:
            desc = desc[77] + "..."
        #return desc + " %d/%d genes hit" % (self.num_hits, len(self.genes))
        return "%s (%d%% of genes show similarity)" % (desc, self.unique_hit_count * 100 / len(self.genes))

    def reverse_strand(self):
        self.reversed = not self.reversed
        self.overall_strand *= -1
        for gene in self.genes:
            gene.reverse()

    def __len__(self):
        return abs(self.start - self.end)

    def _add_label(self, group, v_offset):
        acc = Text(self.accession + '_' + self.ref_cluster_number + ": " + self.full_description, 5, 20 + v_offset)
        if self.accession.startswith('BGC'):
            acc = Text('<a xlink:href="http://mibig.secondarymetabolites.org/repository/' + self.accession + '/index.html#cluster-1" target="_blank">'
                       + self.accession + '</a>: ' + self.full_description, 5, 20 + v_offset)
        elif self.accession.split("_")[0] in get_antismash_db_accessions():
            acc = Text('<a xlink:href="http://antismash-db.secondarymetabolites.org/output/' + self.accession + '/index.html#cluster-' + self.ref_cluster_number[1:] + '" target="_blank">'
                       + self.accession + '_' + self.ref_cluster_number + '</a>: ' + self.full_description, 5, 20 + v_offset)
        acc.set_class("clusterblast-acc")
        group.addElement(acc)
        group.setAttribute('label', self.accession)
        group.setAttribute('description', self.full_description)
        group.set_class('clusterblast-cluster')
        return group

    def get_svg_groups(self, h_offset=0, v_offset=0, scaling=1., screenwidth=1024, colours=None, overview=False):
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
        #Add gene arrows
        arrow_y = line_y + 5
        offset = h_offset - self.start
        if self.reversed:
            offset = h_offset + self.end + self.start
        for i, gene in enumerate(self.genes):
            group = Group()
            arrow = gene.get_arrow_polygon(scaling=scaling, offset=offset,
                               base=arrow_y, colour=colours.get(gene.name, "white"))
            if overview:
                arrow.set_id("%s-%d_all_%d_%d_%s" % ("DUMMY_FILE", self.query_cluster_number, self.query_cluster_number, self.rank, i))
            else:
                arrow.set_id("%s-%d_h%d_%s_%s" % ("DUMMY_FILE", self.query_cluster_number, self.query_cluster_number, self.rank, i))
            group.addElement(arrow)
            #Can be used for domains
            group.set_id("a%s_00%s" % (self.query_cluster_number, i))
            groups.append(group)
        return groups

    @staticmethod
    def from_reference_cluster(cluster, query_cluster_number, score, reference_proteins, rank, num_hits, strand):
        proteins = [reference_proteins[protein] for protein in cluster.proteins]
        cluster = Cluster(query_cluster_number, cluster.cluster_label, cluster.accession, cluster.description, proteins, rank, num_hits, strand)
        for query, subject in score.scored_pairings:
            for gene in cluster.genes:
                if gene.name == subject.name:
                    gene.pairings.append((query, subject))
        return cluster


class QueryCluster(Cluster):
    def __init__(self, cluster_feature):
        cluster_number = cluster_feature.get_cluster_number()
        super().__init__(cluster_number, cluster_number,
                         "%s_%d" % (cluster_feature.parent_record.id, cluster_number),
                         "Query sequence",
                         cluster_feature.cds_children, rank=0)

    def get_svg_groups(self, h_offset=0, v_offset=0, scaling=1., screenwidth=1024, colours=None, overview=False):
        if not colours:
            colours = {}
        groups = []
        group = Group()
        acc = Text(self.description, 5, 20)
        acc.set_class("clusterblast-acc")
        group.addElement(acc)
        line_y = 35 + v_offset
        group.addElement(ShapeBuilder().createLine(10, line_y,
                                                   10 + (screenwidth * 0.75), line_y, strokewidth=1, stroke="grey"))
        group.setAttribute('label', self.description)
        group.set_class('clusterblast-cluster')
        groups.append(group)

        base = line_y + 5
        offset = h_offset + 10 - self.start # 10 for margin
        for n, gene in enumerate(self.genes):
            arrow = gene.get_arrow_polygon(scaling=scaling, offset=offset, base=base, colour=colours.get(gene.name, "white"))
            arrow.set_id("%s-%s_q%s_%s_%s" % ("DUMMY_FILE", self.query_cluster_number, n, self.rank, "all"))
            group.addElement(arrow)
        return groups

class ClusterSVGBuilder:
    def __init__(self, cluster_feature, ranking, reference_proteins):
        if ranking:
            assert reference_proteins
        self.query_cluster = QueryCluster(cluster_feature)
        query_cluster_number = cluster_feature.get_cluster_number()
        self.colour_lookup = build_colour_groups(cluster_feature.cds_children, ranking)
        self.hits = []
        record_prefix = cluster_feature.parent_record.id.split(".", 1)[0]
        num_added = 0
        limit = Config().cb_nclusters
        for cluster, score in ranking:
            if record_prefix == cluster.accession.split("_", 1)[0]:
                continue
            hit_genes = set()
            strand = 0
            for _, subject in score.scored_pairings:
                if subject.name in hit_genes:
                    continue
                hit_genes.add(subject.name)
                if subject.strand == "+":
                    strand += 1
                elif subject.strand == "-":
                    strand -= 1
            if strand < 0:
                strand = -1
            elif strand > 0:
                strand = 1
            cluster = Cluster.from_reference_cluster(cluster, query_cluster_number, score, reference_proteins, num_added + 1, len(hit_genes), strand)
            self.hits.append(cluster)
            num_added += 1
            if num_added >= limit:
                break
        self.max_length = self._size_of_largest_cluster()
        self._organise_strands()

    def get_cluster_descriptions(self):
        return [cluster.full_description for cluster in self.hits]

    def get_cluster_accessions(self):
        return [cluster.accession for cluster in self.hits]

    def _organise_strands(self):
        query_strand = self.query_cluster.overall_strand
        for cluster in self.hits:
            if cluster.overall_strand != query_strand:
                cluster.reverse_strand()

    def _size_of_largest_cluster(self):
        query_length = len(self.query_cluster)
        length = query_length
        for cluster in self.hits:
            if len(cluster) > length:
                length = len(cluster)
        # if this would shrink the query too much, use the minimum allowed
        if query_length / length < Config().cb_min_homology_scale:
            length = query_length / Config().cb_min_homology_scale
        return length

    def get_overview_contents(self, width, height):
        svg = Svg(x=0, y=0, width=width, height=height)
        viewbox = "0 0 %d %d" % (width, height)
        svg.set_viewBox(viewbox)
        svg.set_preserveAspectRatio("none")
        scaling = (width - 20) / self.max_length # -20 for margins
        for group in self.query_cluster.get_svg_groups(h_offset=(self.max_length - len(self.query_cluster)) // 2, scaling=scaling, colours=self.colour_lookup, overview=True):
            svg.addElement(group)
        for n, cluster in enumerate(self.hits):
            for group in cluster.get_svg_groups(v_offset=50 * (n + 1), h_offset=(self.max_length - len(cluster)) // 2,
                             scaling=scaling, colours=self.colour_lookup, overview=True):
                svg.addElement(group)
        return svg.getXML()

    def get_pairing_contents(self, index, width, height):
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
