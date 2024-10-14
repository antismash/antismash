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
from typing import Dict, Iterable, List, Set, Tuple, TypeVar

from antismash.common import secmet

from .data_structures import ReferenceCluster, Score

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


def build_colour_groups(query_cds_features: Iterable[secmet.CDSFeature],
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
