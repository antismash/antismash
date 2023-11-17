# License: GNU Affero General Public License v3 or later
# A copy of GNU AGPL v3 should have been included in this software package in LICENSE.txt.

""" Handles creation of candidate clusters from protoclusters
"""

import bisect
from typing import Iterable, List, Optional

from ...locations import (
    Location,
    connect_locations,
    location_contains_other,
    locations_overlap,
)
from ..protocluster import Protocluster
from .structures import CandidateCluster, CandidateClusterKind


def create_candidates_from_protoclusters(protoclusters: List[Protocluster], circular_wrap_point: int = None,
                                         ) -> List[CandidateCluster]:
    """ Constructs CandidateCluster features from a collection of Protocluster features

        CandidateClusters created may overlap if they are of different kinds.

        Chemical hybrid CandidateClusters may also contain one or more Protoclusters that
        do not share a Protocluster-defining CDS with the others, provided that the
        Protocluster(s) are fully contained within the area covered by other Protoclusters
        that do shared defining CDS features.

        Arguments:
            clusters: a list of Protoluster features
            circular_wrap_point: the coordinate at which a circular record wraps back to the origin,
                                 or None in the case of linear records

        Returns:
            a list of CandidateCluster features
    """
    if not protoclusters:
        return []

    # ensure that only one candidate cluster can be created for each specific combination
    existing: dict[tuple[int, int], CandidateCluster] = {}
    # and allow for keeping some additional singles when handling position conflicts
    singles: set[Protocluster] = set()

    def build_candidates(groups: List[List[Protocluster]], kind: CandidateClusterKind) -> List[CandidateCluster]:
        """ Creates candidates of the given kind, one for each group
        """
        for group in groups:
            assert kind == CandidateClusterKind.SINGLE or len(group) > 1
            candidate = CandidateCluster(kind, sorted(group), circular_wrap_point=circular_wrap_point)
            key = (int(candidate.start), int(candidate.end))
            existing_candidate = existing.get(key)
            if not existing_candidate:
                existing[key] = candidate
                continue

            # handle duplicate locations by checking if any extra protoclusters are in this group
            existing_clusters = set(existing_candidate.protoclusters)
            extras = set(group).difference(existing_clusters)
            if not extras:
                # drop it as a completely redundant, weaker candidate
                continue
            # but if there's extras, promote this group to the existing kind
            replacement = CandidateCluster(existing_candidate.kind, sorted(list(existing_clusters) + list(extras)),
                                           circular_wrap_point=circular_wrap_point)
            # then replace the original
            existing[key] = replacement
            # and track a SINGLE for each of the extras, since they don't quite
            # belong in the stronger candidate grouping
            for extra in extras:
                singles.add(extra)

        return sorted(existing.values())

    unassigned = sorted(protoclusters)  # will contain any protocluster not yet in hybrid/interleaved

    # first all the chemical hybrids
    hybrid_groups, unassigned = _find_hybrids(unassigned, circular_wrap_point)
    candidates = build_candidates(hybrid_groups, CandidateClusterKind.CHEMICAL_HYBRID)

    # then any interleaved
    interleaved_groups, unassigned = _find_interleaved(unassigned, candidates, circular_wrap_point)
    # unassigned won't be updated further, as anything left needs a SINGLE created anyway
    candidates = build_candidates(interleaved_groups, CandidateClusterKind.INTERLEAVED)

    # then neighbouring
    neighbouring_groups = _find_neighbouring(unassigned, candidates)
    candidates = build_candidates(neighbouring_groups, CandidateClusterKind.NEIGHBOURING)

    # finally, create singles for every protocluster not in interleaved or hybrid
    unassigned.extend(singles)  # and those where kind promotion applied
    for cluster in set(unassigned):
        existing_candidate = existing.get((cluster.location.start, cluster.location.end))
        # don't add a single if it's the same coordinates as the parent candidate
        if existing_candidate and cluster in existing_candidate.protoclusters:
            continue
        candidates.append(CandidateCluster(CandidateClusterKind.SINGLE, [cluster],
                                           circular_wrap_point=circular_wrap_point))

    # and as a sanity check, ensure all protoclusters belong to at least one candidate
    assigned = set()
    for candidate in candidates:
        for proto in candidate.protoclusters:
            assigned.add(proto)
    assert len(assigned) == len(protoclusters), f"{len(assigned)} == {len(protoclusters)}"

    return sorted(candidates)


def _merge_sets(groups: Iterable[set[Protocluster]]) -> list[list[Protocluster]]:
    """ Merges all overlapping sets into a larger set. All sets returned are disjoint. """
    groups = [set(group) for group in groups]  # don't modify the inputs
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
    return [sorted(group) for group in ordered if group]  # remove any empty sets


def _find_hybrids(clusters: list[Protocluster], wrap_point: Optional[int],
                  ) -> tuple[list[list[Protocluster]], list[Protocluster]]:
    """ Finds all chemical hybrid groupings within a sorted list of protoclusters """
    unassigned = set(clusters)
    clusters = sorted(clusters, key=lambda x: (x.core_location.start, x.core_location.end))

    # first find all hybrid pairs
    groups = []
    for i, first in enumerate(clusters[:-1]):
        for second in clusters[i+1:]:
            if first.definition_cdses.intersection(second.definition_cdses):
                groups.append({first, second})
                unassigned.discard(first)
                unassigned.discard(second)

    # handle origin-crossing pairs
    first = clusters[0]
    second = clusters[-1]
    if first is not second and first.definition_cdses.intersection(second.definition_cdses):
        groups.append({first, second})
        unassigned.discard(first)
        unassigned.discard(second)

    # merge the pairs into larger groups
    merged_groups = _merge_sets(groups)
    # sort by core_location to ensure vastly different neighbourhood ranges don't cause issues
    clusters = sorted(unassigned, key=lambda x: x.core_location.start)

    def update_if_contained(outer: Location, cluster: Protocluster) -> None:
        if location_contains_other(outer, cluster.core_location):
            group.append(cluster)
            unassigned.discard(cluster)

    # then find any unassigned cluster that is fully contained by a hybrid group
    for group in merged_groups:
        core = connect_locations([proto.core_location for proto in group], wrap_point=wrap_point)
        start = bisect.bisect_left([cluster.core_start for cluster in clusters], core.start) - 1
        index = max(0, start)
        for cluster in clusters[index:]:
            if cluster.location.start > core.end:
                break
            update_if_contained(core, cluster)
        if len(core.parts) > 1:
            for cluster in clusters:
                if cluster.location.start > core.parts[-1].end:
                    break
                update_if_contained(core, cluster)
    return [sorted(group) for group in merged_groups], sorted(unassigned)


def _find_interleaved_candidates(candidates: list[CandidateCluster]) -> list[set[Protocluster]]:
    """ Returns disjoint sets of candidates with any protoclusters that overlap core locations """
    groups = []

    # start with any interleaved candidates
    for i, candidate in enumerate(candidates[:-1]):
        for other_candidate in candidates[i+1:]:
            if locations_overlap(candidate.core_location, other_candidate.core_location):
                groups.append(set(candidate.protoclusters + other_candidate.protoclusters))

    # handle origin-crossing pairs
    if len(candidates) > 1:
        candidate = candidates[0]
        other_candidate = candidates[-1]
        if locations_overlap(candidate.core_location, other_candidate.core_location):
            groups.append(set(candidate.protoclusters + other_candidate.protoclusters))
    return groups


def _find_cross_origin_interleaved(candidates: list[CandidateCluster], unassigned: list[Protocluster],
                                   existing_groups: list[set[Protocluster]], wrap_point: Optional[int],
                                   ) -> set[Protocluster]:
    """ Returns a set of protoclusters with cores that overlap with the given candidates,
        adding any overlaps to the existing groups.
    """
    found: set[Protocluster] = set()
    if not (unassigned and candidates):
        return found

    # it's possible that only the non-core areas cross the origin, in which case
    # there's nothing to check
    if not any(candidate.core_crosses_origin() for candidate in candidates):
        return found

    core = connect_locations([c.core_location for c in candidates if c.core_crosses_origin()], wrap_point)
    core_group: set[Protocluster] = set()
    for candidate in candidates:
        if candidate.core_crosses_origin():
            core_group.update([proto for proto in candidate.protoclusters if proto.core_location.crosses_origin()])
    assert core_group
    for direction in [-1, 1]:
        index = 0 if direction == 1 else -1  # don't check 0 twice, especially if it ends the loop early
        while (abs(index) < len(unassigned)  # don't go out of bounds
               and len(found) < len(unassigned)):  # avoid cases of whole-record overlaps being infinite
            cluster = unassigned[index]
            # once disconnected candidates are hit, no more interleaved are possible
            if not locations_overlap(cluster.core_location, core):
                break
            core_group.add(cluster)
            found.add(cluster)
            index += direction
    if any(core_group == set(candidate.protoclusters) for candidate in candidates):
        return set()
    if len(core_group) > 1:
        existing_groups.append(core_group)
    return found


def _find_interleaved(clusters: list[Protocluster], candidates: list[CandidateCluster],
                      wrap_point: Optional[int]) -> tuple[list[list[Protocluster]], list[Protocluster]]:
    """ Finds all interleaved groups from combinations of existing candidates
        and single protoclusters (both inputs must be sorted)
    """
    found = set()
    groups = _find_interleaved_candidates(candidates)

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

    # and origin-crossing combinations of unassigned and candidates
    found.update(_find_cross_origin_interleaved(candidates, unassigned_by_core, groups, wrap_point=wrap_point))

    return _merge_sets(groups), sorted(set(clusters).difference(found))


def _find_neighbouring_candidates(candidates: list[CandidateCluster]) -> list[set[Protocluster]]:
    """ Returns sets of candidates with overlapping locations """
    groups = []
    # candidates overlapping with candidates
    for i, first_candidate in enumerate(candidates[:-1]):
        for second_candidate in candidates[i+1:]:
            if not first_candidate.overlaps_with(second_candidate):
                continue
            groups.append(set(first_candidate.protoclusters).union(set(second_candidate.protoclusters)))

    # handle origin-crossing pairs
    if len(candidates) > 1:
        candidate = candidates[0]
        i = -1
        while 0 < i < len(candidates):
            other_candidate = candidates[i]
            if locations_overlap(candidate.location, other_candidate.location):
                groups.append(set(candidate.protoclusters + other_candidate.protoclusters))
            i += 1
    return groups


def _find_neighbouring_protoclusters(protoclusters: list[Protocluster]) -> list[set[Protocluster]]:
    """ Returns sets of protoclusters with overlapping locations """
    groups = []
    for i, first_cluster in enumerate(protoclusters[:-1]):
        for second_cluster in protoclusters[i+1:]:
            if not first_cluster.overlaps_with(second_cluster):
                continue
            groups.append({first_cluster, second_cluster})

    # again, origin-crossing needs to be handled between protoclusters
    if len(protoclusters) > 1:
        first_cluster = protoclusters[0]
        second_cluster = protoclusters[-1]
        if first_cluster is not second_cluster and first_cluster.overlaps_with(second_cluster):
            groups.append({first_cluster, second_cluster})
    return groups


def _find_neighbouring(singles: list[Protocluster], candidates: list[CandidateCluster]) -> list[list[Protocluster]]:
    """ Finds all neighbouring groups from combinations of existing candidates
        and single protoclusters (both inputs must be sorted)
    """
    groups = _find_neighbouring_candidates(candidates)
    unassigned = set(singles)

    # singles overlapping with candidates
    for single in singles:
        index = max(0, bisect.bisect_left(candidates, single) - 1)
        for candidate in (candidates[index:] + candidates[:1]):
            if candidate.location.start > single.location.end:
                break
            if single.overlaps_with(candidate):
                groups.append(set(candidate.protoclusters).union({single}))
                unassigned.discard(single)
    # and origin-crossing combinations of singles and candidates
    edges = []
    if unassigned and candidates:
        if candidates[0].crosses_origin():
            edges.append(candidates[0])
        if len(candidates) > 1 and candidates[-1].crosses_origin():
            edges.append(candidates[-1])

        for candidate in edges:
            for single in unassigned:
                if locations_overlap(single.location, candidate.location):
                    groups.append(set(candidate.protoclusters + (single,)))
                    break
    groups.extend(_find_neighbouring_protoclusters(sorted(unassigned)))
    return _merge_sets(groups)
