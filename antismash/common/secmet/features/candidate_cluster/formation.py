# License: GNU Affero General Public License v3 or later
# A copy of GNU AGPL v3 should have been included in this software package in LICENSE.txt.

""" Handles creation of candidate clusters from protoclusters
"""

import bisect
from typing import Iterable, List, Set, Tuple

from ...locations import FeatureLocation, location_contains_other, locations_overlap
from ..protocluster import Protocluster
from .structures import CandidateCluster, CandidateClusterKind


def create_candidates_from_protoclusters(protoclusters: List[Protocluster]) -> List[CandidateCluster]:
    """ Constructs CandidateCluster features from a collection of Protocluster features

        CandidateClusters created may overlap if they are of different kinds.

        Chemical hybrid CandidateClusters may also contain one or more Protoclusters that
        do not share a Protocluster-defining CDS with the others, provided that the
        Protocluster(s) are fully contained within the area covered by other Protoclusters
        that do shared defining CDS features.

        Arguments:
            clusters: a list of Protoluster features

        Returns:
            a list of CandidateCluster features
    """
    if not protoclusters:
        return []

    # ensure that only one candidate cluster can be created for each specific combination
    existing: set[tuple[int, int, tuple[str, ...]]] = set()

    def build_candidates(groups: List[List[Protocluster]], kind: CandidateClusterKind) -> List[CandidateCluster]:
        """ Creates candidates of the given kind, one for each group
        """
        candidates = []
        for group in groups:
            assert kind == CandidateClusterKind.SINGLE or len(group) > 1
            candidate = CandidateCluster(kind, sorted(group))
            key = (candidate.location.start, candidate.location.end, tuple(sorted(candidate.products)))
            if key not in existing:
                candidates.append(candidate)
                existing.add(key)
        return candidates

    unassigned = sorted(protoclusters)  # will contain any protocluster not yet in hybrid/interleaved

    # first all the chemical hybrids
    hybrid_groups, unassigned = _find_hybrids(unassigned)
    candidates = build_candidates(hybrid_groups, CandidateClusterKind.CHEMICAL_HYBRID)

    # then any interleaved
    interleaved_groups, unassigned = _find_interleaved(unassigned, candidates)
    # unassigned won't be updated further, as anything left needs a SINGLE created anyway
    new = build_candidates(interleaved_groups, CandidateClusterKind.INTERLEAVED)
    candidates = sorted(candidates + new)  # since find_neighbouring requires sorted by location

    # then neighbouring
    neighbouring_groups = _find_neighbouring(unassigned, candidates)
    new = build_candidates(neighbouring_groups, CandidateClusterKind.NEIGHBOURING)
    candidates.extend(new)  # sorting postponed, since it's not critical for next step

    # finally, create singles for every protocluster not in interleaved or hybrid
    for cluster in unassigned:
        if (cluster.location.start, cluster.location.end) not in existing:
            candidates.append(CandidateCluster(CandidateClusterKind.SINGLE, [cluster]))

    # and as a sanity check, ensure all protoclusters belong to at least one candidate
    assigned = set()
    for candidate in candidates:
        for proto in candidate.protoclusters:
            assigned.add(proto)
    assert len(assigned) == len(protoclusters), f"{len(assigned)} == {len(protoclusters)}"

    return sorted(candidates)  # again, reorder by location


def _merge_sets(groups: Iterable[Set[Protocluster]]) -> List[List[Protocluster]]:
    """ Merges all overlapping sets into a larger set. All sets returned are disjoint. """
    # once ordered by earliest start location, a single pass should be sufficient
    ordered = sorted(groups, key=lambda group: min(cluster.location.start for cluster in group))
    for i, first in enumerate(ordered[:-1]):
        if not first:
            continue
        for second in ordered[i+1:]:
            if first.isdisjoint(second):
                continue
            # merge into the earlier group and make the second unable to hit any other
            first.update(second)
            second.clear()
    return [sorted(group) for group in ordered if group]


def _find_hybrids(clusters: List[Protocluster]) -> Tuple[List[List[Protocluster]], List[Protocluster]]:
    """ Finds all chemical hybrid groupings within a sorted list of protoclusters """
    unassigned = set(clusters)
    clusters = sorted(clusters, key=lambda x: (x.core_location.start, x.core_location.end))

    # first find all hybrid pairs
    groups = []
    for i, first in enumerate(clusters[:-1]):
        for second in clusters[i+1:]:
            if not first.overlaps_with(second):
                break
            if first.definition_cdses.intersection(second.definition_cdses):
                groups.append({first, second})
                unassigned.discard(first)
                unassigned.discard(second)

    # merge the pairs into larger groups
    merged_groups = _merge_sets(groups)
    # sort by core_location to ensure vastly different neighbourhood ranges don't cause issues
    clusters = sorted(unassigned, key=lambda x: x.core_location.start)

    # then find any unassigned cluster that is fully contained by a hybrid group
    for group in merged_groups:
        core = FeatureLocation(min(proto.core_location.start for proto in group),
                               max(proto.core_location.end for proto in group))
        index = max(0, bisect.bisect_left([cluster.core_location.start for cluster in clusters], core.start) - 1)
        for cluster in clusters[index:]:
            if cluster.location.start > core.end:
                break
            if location_contains_other(core, cluster.core_location):
                group.append(cluster)
                unassigned.discard(cluster)

    return [sorted(group) for group in merged_groups], sorted(unassigned)


def _find_interleaved(clusters: List[Protocluster], candidates: List[CandidateCluster]
                      ) -> Tuple[List[List[Protocluster]], List[Protocluster]]:
    """ Finds all interleaved groups from combinations of existing candidates
        and single protoclusters (both inputs must be sorted)
    """
    found = set()
    groups = []

    # start with any interleaved candidates
    for i, candidate in enumerate(candidates[:-1]):
        for other_candidate in candidates[i+1:]:
            if locations_overlap(candidate.core_location, other_candidate.core_location):
                groups.append(set(candidate.protoclusters + other_candidate.protoclusters))

    # sort unassigned by core location, as that's the important part
    unassigned_by_core = sorted(clusters, key=lambda x: x.core_location.start)

    # unassigned overlapping with unassigned
    for i, cluster in enumerate(unassigned_by_core):
        for other_cluster in unassigned_by_core[i + 1:]:
            if cluster.core_location.end <= other_cluster.core_location.start:
                break
            if locations_overlap(cluster.core_location, other_cluster.core_location):
                groups.append({cluster, other_cluster})
                found.add(cluster)
                found.add(other_cluster)

    # unassigned overlapping with candidates
    for cluster in unassigned_by_core:
        index = max(0, bisect.bisect_left(candidates, cluster) - 1)
        for candidate in candidates[index:]:
            if candidate.location.start > cluster.location.end:
                break
            if locations_overlap(candidate.core_location, cluster.core_location):
                groups.append(set(candidate.protoclusters + (cluster,)))
                found.add(cluster)

    return _merge_sets(groups), sorted(set(clusters).difference(found))


def _find_neighbouring(singles: List[Protocluster], candidates: List[CandidateCluster]) -> List[List[Protocluster]]:
    """ Finds all neighbouring groups from combinations of existing candidates
        and single protoclusters (both inputs must be sorted)
    """
    groups = []
    # candidates overlapping with candidates
    for i, first_candidate in enumerate(candidates[:-1]):
        for second_candidate in candidates[i+1:]:
            if not first_candidate.overlaps_with(second_candidate):
                break
            groups.append(set(first_candidate.protoclusters).union(set(second_candidate.protoclusters)))

    # singles overlapping with candidates
    for single in singles:
        index = max(0, bisect.bisect_left(candidates, single) - 1)
        for candidate in candidates[index:]:
            if candidate.location.start > single.location.end:
                break
            if single.overlaps_with(candidate):
                groups.append(set(candidate.protoclusters).union({single}))

    # singles overlapping with singles
    for i, first_single in enumerate(singles[:-1]):
        for second_single in singles[i+1:]:
            if not first_single.overlaps_with(second_single):
                break
            groups.append({first_single, second_single})

    return _merge_sets(groups)
