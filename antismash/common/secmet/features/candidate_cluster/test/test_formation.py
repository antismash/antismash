# License: GNU Affero General Public License v3 or later
# A copy of GNU AGPL v3 should have been included in this software package in LICENSE.txt.

# for test files, silence irrelevant and noisy pylint warnings
# pylint: disable=use-implicit-booleaness-not-comparison,protected-access,missing-docstring

from functools import partial
import unittest

from antismash.common.secmet.features import CandidateCluster, Protocluster
from antismash.common.secmet.features.candidate_cluster import (
    create_candidates_from_protoclusters as creator,
)
from antismash.common.secmet.features.candidate_cluster.formation import (
    _find_cross_origin_interleaved as find_origin_interleaved,
    _find_hybrids as find_hybrids,
    _find_interleaved as find_interleaved,
    _find_neighbouring as find_neighbouring,
    _merge_sets as merge,
)
from antismash.common.secmet.locations import (
    CompoundLocation,
    FeatureLocation,
    location_bridges_origin,
)
from antismash.common.secmet.qualifiers.gene_functions import GeneFunction
from antismash.common.secmet.test.helpers import DummyCandidateCluster, DummyCDS

from .test_candidate_cluster import create_cds, create_cluster


def create_compound_loc(pairs: list[tuple[int, int]]) -> CompoundLocation:
    assert len(pairs) > 1
    return CompoundLocation([FeatureLocation(start, end, 1) for start, end in pairs])


def make_core_cds(location, products):
    assert products
    cds = DummyCDS(location=location)
    for product in products:
        cds.gene_functions.add(GeneFunction.CORE, tool="tool", description="dummy", product=product)
    return cds


class TestCreation(unittest.TestCase):
    def test_creation_empty(self):
        assert creator([]) == []

    def test_creation_single(self):
        created = creator([create_cluster(3, 8, 71, 76, 'a')])
        print(created)
        assert len(created) == 1
        assert created[0].location == FeatureLocation(3, 76, 1)
        assert created[0].kind == CandidateCluster.kinds.SINGLE

    def test_creation_neighbours(self):
        cluster = create_cluster(3, 8, 71, 76, 'a')
        extra_cluster = create_cluster(50, 100, 120, 170, 'b')
        created = creator([cluster, extra_cluster])
        print(created)
        assert len(created) == 3
        expected_location = FeatureLocation(cluster.location.start, extra_cluster.location.end, 1)
        assert created[0].kind == CandidateCluster.kinds.NEIGHBOURING and created[0].location == expected_location
        assert created[1].kind == CandidateCluster.kinds.SINGLE and created[1].location == cluster.location
        assert created[2].kind == CandidateCluster.kinds.SINGLE and created[2].location == extra_cluster.location

    def test_creation_coreoverlap(self):
        cluster = create_cluster(3, 8, 71, 76, 'a')
        extra_cluster = create_cluster(50, 60, 120, 170, 'b')
        # create a CDS within both clusters that has a product from only one cluster
        cds = create_cds(60, 65, ["a"])
        cluster.add_cds(cds)
        extra_cluster.add_cds(cds)

        created = creator([cluster, extra_cluster])
        print(created)
        assert len(created) == 1
        candidate_cluster = created[0]
        assert candidate_cluster.kind == CandidateCluster.kinds.INTERLEAVED
        assert candidate_cluster.location == FeatureLocation(3, 170, 1)

    def test_creation_hybrid(self):
        cluster = create_cluster(3, 8, 71, 76, 'a')
        hybrid_cluster = create_cluster(50, 60, 120, 170, 'b')

        # insert the cds that will cause the hybrid call
        cds_ab = create_cds(60, 65, ["a", "b"])
        cluster.add_cds(cds_ab)
        hybrid_cluster.add_cds(cds_ab)

        created = creator([cluster, hybrid_cluster])
        print(created)
        assert len(created) == 1
        candidate_cluster = created[0]
        assert candidate_cluster.kind == CandidateCluster.kinds.CHEMICAL_HYBRID
        assert candidate_cluster.location == FeatureLocation(3, 170, 1)

    def test_creation_mixed(self):
        cluster = create_cluster(3, 8, 71, 76, 'a')
        hybrid_cluster = create_cluster(50, 60, 120, 170, 'b')
        overlap_cluster = create_cluster(80, 90, 130, 180, 'o')
        neighbour_cluster = create_cluster(50, 210, 260, 270, 'a')
        isolated_cluster = create_cluster(450, 500, 550, 600, 'alone')

        # insert the cds that will cause the hybrid call
        cds_ab = create_cds(60, 65, ["a", "b"])
        cluster.add_cds(cds_ab)
        hybrid_cluster.add_cds(cds_ab)

        created = creator([cluster, hybrid_cluster, overlap_cluster, neighbour_cluster, isolated_cluster])
        print(created)

        assert len(created) == 5
        assert created[0].location == FeatureLocation(3, 270, 1)
        assert created[0].kind == CandidateCluster.kinds.NEIGHBOURING
        assert created[0].protoclusters == (cluster, neighbour_cluster, hybrid_cluster, overlap_cluster)

        assert created[1].location == FeatureLocation(3, 180, 1)
        assert created[1].kind == CandidateCluster.kinds.INTERLEAVED
        assert created[1].protoclusters == (cluster, hybrid_cluster, overlap_cluster)

        assert created[2].location == FeatureLocation(3, 170, 1)
        assert created[2].kind == CandidateCluster.kinds.CHEMICAL_HYBRID
        assert created[2].protoclusters == (cluster, hybrid_cluster)

        assert created[3].location == FeatureLocation(50, 270, 1)
        assert created[3].kind == CandidateCluster.kinds.SINGLE
        assert created[3].protoclusters == (neighbour_cluster,)

        assert created[4].location == FeatureLocation(450, 600, 1)
        assert created[4].kind == CandidateCluster.kinds.SINGLE
        assert created[4].protoclusters == (isolated_cluster,)

    def test_hybrid_interactions(self):
        cluster = create_cluster(3, 8, 171, 176, "a")
        hybrid = create_cluster(3, 8, 50, 55, "b")
        contained = create_cluster(80, 90, 100, 110, "c")  # will form part of hybrid

        hybrid_cds = create_cds(8, 50, ["a", "b"])
        cluster.add_cds(hybrid_cds)
        hybrid.add_cds(hybrid_cds)

        for overlapping in [create_cluster(120, 130, 200, 250, "d"),
                            create_cluster(60, 70, 200, 250, "d")]:

            created = creator([cluster, hybrid, contained, overlapping])

            assert len(created) == 2
            assert created[0].location == FeatureLocation(3, 250, 1)
            assert created[0].kind == CandidateCluster.kinds.INTERLEAVED
            assert created[0].protoclusters == tuple(sorted([cluster, hybrid, contained, overlapping]))

            assert created[1].location == FeatureLocation(3, 176, 1)
            assert created[1].kind == CandidateCluster.kinds.CHEMICAL_HYBRID
            assert created[1].protoclusters == (cluster, hybrid, contained)

    def test_interleaving(self):
        # these first two hybrid clumps should be interleaved
        first_hybrid_clusters = [create_cluster(30, 60, 120, 150, "a"),
                                 create_cluster(60, 90, 150, 180, "b")]
        cds = create_cds(90, 120, ["a", "b"])
        for cluster in first_hybrid_clusters:
            cluster.add_cds(cds)

        second_hybrid_clusters = [create_cluster(90, 120, 250, 280, "c"),
                                  create_cluster(190, 220, 280, 310, "d")]
        cds = create_cds(220, 250, ["c", "d"])
        for cluster in second_hybrid_clusters:
            cluster.add_cds(cds)

        # this non-hybrid should also be included in the interleaved
        single = create_cluster(230, 250, 410, 430, "e")
        # this hybrid should not
        standalone = [create_cluster(1000, 1100, 1400, 1500, "f"),
                      create_cluster(1100, 1200, 1500, 1600, "g")]
        cds = create_cds(1300, 1400, ["f", "g"])
        for cluster in standalone:
            cluster.add_cds(cds)

        created = creator(first_hybrid_clusters + second_hybrid_clusters + [single] + standalone)

        assert len(created) == 4
        assert created[0].location == FeatureLocation(30, 430, 1)
        assert created[0].core_location == FeatureLocation(60, 410, 1)
        assert created[0].kind == CandidateCluster.kinds.INTERLEAVED
        assert created[0].protoclusters == tuple(first_hybrid_clusters + second_hybrid_clusters + [single])

        assert created[1].location == FeatureLocation(30, 180, 1)
        assert created[1].protoclusters == tuple(first_hybrid_clusters)
        assert created[2].location == FeatureLocation(90, 310, 1)
        assert created[2].protoclusters == tuple(second_hybrid_clusters)
        for cand in created[1:3]:
            assert cand.kind == CandidateCluster.kinds.CHEMICAL_HYBRID

        assert created[3].location == FeatureLocation(1000, 1600, 1)
        assert created[3].kind == CandidateCluster.kinds.CHEMICAL_HYBRID

    def test_interleaving_order(self):
        clusters = [create_cluster(1000, 1100, 1400, 1500, "a"),
                    create_cluster(1050, 2000, 3000, 4000, "b"),  # sorts second due to neighbouring
                    create_cluster(1100, 1200, 1500, 1600, "c")]
        assert sorted(clusters) == clusters
        created = creator(clusters)
        assert len(created) == 3
        assert created[0].kind == CandidateCluster.kinds.NEIGHBOURING
        assert created[0].location == FeatureLocation(1000, 4000, 1)
        assert created[1].kind == CandidateCluster.kinds.INTERLEAVED
        assert created[1].location == FeatureLocation(1000, 1600, 1)
        assert created[2].kind == CandidateCluster.kinds.SINGLE
        assert created[2].location == FeatureLocation(1050, 4000, 1)

    def test_contained_neighbours(self):
        # if a larger protocluster contains another entirely in its neighbourhood
        # it shouldn't create a single which has the exact same coordinates
        # as the neighbour
        big = create_cluster(0, 100, 130, 230, "a")
        small = create_cluster(10, 20, 30, 40, "b")

        created = creator([big, small])

        assert len(created) == 2
        assert created[0].kind == CandidateCluster.kinds.NEIGHBOURING
        assert created[0].location == big.location
        assert created[0].products == ["a", "b"]

        assert created[1].kind == CandidateCluster.kinds.SINGLE
        assert created[1].location == small.location
        assert created[1].products == ["b"]

    def test_multiple_contained_in_hybrid(self):
        clusters = [
            create_cluster(142916, 162916, 191227, 209850, "a"),  # covers whole region
            create_cluster(142916, 162916, 169708, 189708, "b"),  # causes chem hybrid
            create_cluster(154261, 174261, 175881, 195881, "c"),  # interleaved
            create_cluster(160464, 180464, 183155, 203155, "d"),  # interleaved
        ]
        cds = create_cds(162916, 169708, ["a", "b"])
        clusters[0].add_cds(cds)
        clusters[1].add_cds(cds)

        created = creator(clusters)
        assert len(created) == 1
        assert created[0].kind == CandidateCluster.kinds.CHEMICAL_HYBRID
        assert created[0].protoclusters == tuple(clusters)

    def test_overlap_interleave(self):
        # mimics an edge case where an interleaved was formed due to a small overlap
        # in two genes, and the interleaved was dropped because it had the same
        # coordinates, despite having different products
        clusters = [
            create_cluster(100, 500, 800, 1000, "hybA"),
            create_cluster(100, 500, 800, 1000, "hybB"),
            create_cluster(700, 790, 900, 1000, "inter"),
        ]
        cds = create_cds(500, 800, ["hybA", "hybB"])
        clusters[0].add_cds(cds)
        clusters[1].add_cds(cds)
        cds = create_cds(790, 900, ["inter"])
        clusters[2].add_cds(cds)
        created = creator(clusters)

        assert created[0].kind == CandidateCluster.kinds.CHEMICAL_HYBRID
        assert sorted(created[0].protoclusters, key=id) == sorted(clusters, key=id)

        assert created[1].kind == CandidateCluster.kinds.SINGLE
        assert created[1].protoclusters == (clusters[2],)

        assert len(created) == 2

    def test_protocluster_promotion(self):
        # mimics a case where a neighbouring/interleaved protocluster was lost
        # due to a small overlap in core genes, while near a contig edge
        # and the other clusters formed a hybrid
        contig_edge = 2300
        clusters = [
            create_cluster(0, 200, 800, 1000, "a"),
            create_cluster(500, 750, 1500, contig_edge, "b"),
            create_cluster(1750, 2100, 2200, contig_edge, "c"),  # would be neighbouring
        ]
        # defining genes
        hybrid = create_cds(750, 800, ["a", "b"])
        clusters[0].add_cds(hybrid)
        clusters[1].add_cds(hybrid)

        created = creator(clusters)
        assert len(created) == 2, f"extras {created[2:]}"
        assert created[0].kind == CandidateCluster.kinds.CHEMICAL_HYBRID
        assert created[0].products == ["a", "b", "c"]
        assert created[1].kind == CandidateCluster.kinds.SINGLE
        assert created[1].products == ["c"]

    def test_real_over_origin(self):
        # a quite complicated setup, with neighbourhoods and cores crossing the origin
        # complete with a chemical hybrid in the cross-origin cores
        #   --NO-|-
        #    ---C|OR1---
        #    ---C|ORE2---
        #        |   -UN-
        neighbour_over = create_cluster(0, 1, 2, 3, "neighbour_over")
        neighbour_over.location = create_compound_loc([(59, 100), (0, 6)])
        neighbour_over.core_location = FeatureLocation(79, 80, 1)

        short_core_over = create_cluster(0, 1, 2, 3, "short_core_over")
        short_core_over.location = create_compound_loc([(70, 100), (0, 28)])
        short_core_over.core_location = create_compound_loc([(90, 100), (0, 8)])

        long_core_over = create_cluster(0, 1, 2, 3, "long_core_over")
        long_core_over.location = create_compound_loc([(70, 100), (0, 39)])
        long_core_over.core_location = create_compound_loc([(90, 100), (0, 19)])

        unbroken = create_cluster(16, 26, 28, 38, "unbroken")

        # the two cross-origin cores with the same start need to share that start cds
        shared_core = DummyCDS()
        short_core_over._definition_cdses = {shared_core}
        long_core_over._definition_cdses = {shared_core}

        created = creator([neighbour_over, short_core_over, long_core_over, unbroken], circular_wrap_point=100)
        assert len(created) == 4
        hybrids = [c for c in created if c.kind == c.kinds.CHEMICAL_HYBRID]
        assert len(hybrids) == 1, hybrids
        hybrid = hybrids[0]
        assert set(hybrid.products) == {short_core_over.product, long_core_over.product}
        assert hybrid.start == short_core_over.start == long_core_over.start
        assert hybrid.end == long_core_over.end

        interleaved = [c for c in created if c.kind == c.kinds.INTERLEAVED]
        assert not interleaved

        neighbours = [c for c in created if c.kind == c.kinds.NEIGHBOURING]
        assert len(neighbours) == 1
        neighbour = neighbours[0]
        assert set(neighbour.products) == {neighbour_over.product, unbroken.product,
                                           short_core_over.product, long_core_over.product}
        assert neighbour.start == neighbour_over.start
        assert neighbour.end == long_core_over.end

        singles = [c for c in created if c.kind == c.kinds.SINGLE]
        assert len(singles) == 2
        assert {c.protoclusters[0].product for c in singles} == {neighbour_over.product, unbroken.product}

    def test_real_whole_contig(self):
        # another complicated case, found in NZ_AP017310.1, where a circular plasmid is almost
        # entirely covered by four protoclusters, with three of them being over the origin
        #
        # the candidate formed must also be over the origin instead of (0, record length)

        record_length = 98206

        cdses = [
            make_core_cds(FeatureLocation(85398, 89805, -1), ["a", "b"]),
            make_core_cds(CompoundLocation([
                FeatureLocation(0, 3327, -1),
                FeatureLocation(95431, record_length, -1),
            ]), ["a", "b"]),
            make_core_cds(FeatureLocation(15291, 22440, -1), ["a", "b", "c"]),
            make_core_cds(FeatureLocation(6497, 15299, -1), ["a", "b"]),
            make_core_cds(FeatureLocation(27053, 32159, -1), ["single"]),
        ]

        shared_core = CompoundLocation([
            FeatureLocation(85398, record_length, 1),
            FeatureLocation(0, 22440, 1),
        ])
        shared_neighbourhood = CompoundLocation([
            FeatureLocation(65398, record_length, 1),
            FeatureLocation(0, 42440, 1),
        ])
        first = Protocluster(shared_core, shared_neighbourhood, "tool", "a", 20, 20, "a")
        second = Protocluster(shared_core, shared_neighbourhood, "tool", "b", 20, 20, "b")
        third = Protocluster(
            FeatureLocation(15291, 32000, -1),
            CompoundLocation([FeatureLocation(93497, record_length, 1), FeatureLocation(0, 52159, 1)]),
            "tool", "c", 20, 20, "c",
        )
        non_hybrid = create_cluster(7053, 27053, 32159, 52159, "single")
        for proto in [first, second, third, non_hybrid]:
            for cds in cdses:
                if cds.is_contained_by(proto):
                    proto.add_cds(cds)
            assert proto.definition_cdses

        found = creator([first, second, third, non_hybrid], circular_wrap_point=record_length)
        assert len(found) == 2
        chemical, single = found

        assert chemical.kind == CandidateCluster.kinds.CHEMICAL_HYBRID
        assert chemical.crosses_origin()
        assert chemical.location == CompoundLocation([
            FeatureLocation(65398, record_length, 1),
            FeatureLocation(0, 52159, 1),
        ])

        # based only on cores, there would be an interleaved,
        # but the interleaved candidate would have the same coordinates as the hybrid
        # and is discarded

        assert single.kind == CandidateCluster.kinds.SINGLE
        assert not single.crosses_origin()
        assert single.location == non_hybrid.location


class TestHelpers(unittest.TestCase):
    def setUp(self):
        self.record_length = 1000
        self.first = create_cluster(800, 850, 900, 950, "first")
        core = CompoundLocation([FeatureLocation(850, 1000, 1), FeatureLocation(0, 100, 1)])
        surrounds = CompoundLocation([FeatureLocation(800, 1000, 1), FeatureLocation(0, 200, 1)])
        self.core_over_origin = Protocluster(core, surrounds, "tool", "core_over_origin", 50, 100, "rule")
        self.isolated = create_cluster(400, 500, 600, 700, "isolated")
        self.unassigned = [self.first, self.core_over_origin, self.isolated]

    def test_hybrids(self):
        shared_cds = DummyCDS(850, 900)

        # no groups will form because definition CDSes aren't shared, despite overlaps
        groups, leftovers = find_hybrids([self.first, self.core_over_origin], self.record_length)
        assert not groups
        assert leftovers

        self.first._definition_cdses.add(shared_cds)
        assert self.first.definition_cdses == {shared_cds}
        self.core_over_origin._definition_cdses.add(shared_cds)
        assert self.first.definition_cdses == self.core_over_origin.definition_cdses

        assert self.core_over_origin.crosses_origin()
        # now that 1 and 2 share a core gene, a group will form
        groups, leftovers = find_hybrids(self.unassigned, self.record_length)
        assert len(groups) == 1
        assert set(groups[0]) == {self.first, self.core_over_origin}
        assert leftovers == [self.isolated]

    def test_neighbouring(self):
        neighbours = find_neighbouring(self.unassigned, [])
        assert len(neighbours) == 1
        assert set(neighbours[0]) == {self.first, self.core_over_origin}

    def test_neighbouring_candidates(self):
        extra = create_cluster(700, 750, 800, 850, "a")
        candidates = [DummyCandidateCluster(clusters=[extra])]
        neighbours = find_neighbouring(self.unassigned, candidates)
        assert len(neighbours) == 1
        assert set(neighbours[0]) == {extra, self.first, self.core_over_origin}

    def test_neighbouring_cross_origin_candidates(self):
        extra = create_cluster(700, 750, 800, 850, "extra")
        neighbourhood_over_origin = create_cluster(0, 20, 50, 80, "extra_over")
        # override the neighbourhood to extend over the origin
        neighbourhood_over_origin.location = CompoundLocation([
            FeatureLocation(980, 1000, 1),
            FeatureLocation(0, 80, 1)
        ])
        candidates = [
            DummyCandidateCluster(clusters=[self.core_over_origin]),
            DummyCandidateCluster(clusters=[neighbourhood_over_origin])
        ]
        neighbours = find_neighbouring([self.first, self.isolated, extra], candidates)
        assert len(neighbours) == 1
        assert set(neighbours[0]) == {extra, self.first, self.core_over_origin, neighbourhood_over_origin}

    def test_find_interleaved(self):
        groups, leftovers = find_interleaved(self.unassigned, [], self.record_length)
        assert len(groups) == 1
        assert set(groups[0]) == {self.first, self.core_over_origin}
        assert leftovers == [self.isolated]

    def test_cross_origin_interleaved(self):
        core = FeatureLocation(50, 100)
        surrounds = CompoundLocation([FeatureLocation(950, 1000, 1), FeatureLocation(0, 150, 1)])
        surrounds_over_origin = Protocluster(core, surrounds, "tool", "surrounds", 50, 100, "rule")
        candidates = sorted([
            DummyCandidateCluster(clusters=[surrounds_over_origin, self.core_over_origin]),
            DummyCandidateCluster(clusters=[surrounds_over_origin]),
        ])
        assert candidates[0].crosses_origin() and location_bridges_origin(candidates[0].core_location)
        assert candidates[1].crosses_origin() and not location_bridges_origin(candidates[1].core_location)
        assert all(candidate.crosses_origin() for candidate in candidates)
        existing = [{surrounds_over_origin}, {self.core_over_origin}]  # this will be modified by the function
        original = [set(group) for group in existing]
        found = find_origin_interleaved(candidates, sorted([self.first, self.isolated]), existing, self.record_length)
        assert len(found) == 1
        assert found == {self.first}
        assert existing != original
        assert existing == [{surrounds_over_origin}, {self.core_over_origin},
                            {self.core_over_origin, self.first}]

    def test_cross_origin_real(self):
        # NZ_CP029364.1 version 28-MAR-2024
        # crashed with two cross-origin chemical hybrid clusters and an unrelated, non-overlapping cluster
        # one of the cross-origin protoclusters was left alone in a group after hybrids were formed,
        # despite being in a chemical hybrid and not near the independent cluster
        first = Protocluster(
            CompoundLocation([FeatureLocation(4124436, 4154245, 1), FeatureLocation(0, 7971, 1)]),
            CompoundLocation([FeatureLocation(4104436, 4154245, 1), FeatureLocation(0, 27971, 1)]),
            "tool", "NRPS", 20, 20, "NRPS",
        )
        second = Protocluster(
            FeatureLocation(4146, 16139, 1),
            CompoundLocation([FeatureLocation(4148391, 4154245, 1), FeatureLocation(0, 26139, 1)]),
            "tool", "betalactone", 10, 10, "betalactone",
        )

        shared = make_core_cds(FeatureLocation(4146, 7971, strand=1), ["NRPS", "betalactone"])
        first.add_cds(shared)
        second.add_cds(shared)

        assert set(first.definition_cdses).intersection(second.definition_cdses)

        third = Protocluster(
            FeatureLocation(119498, 214508, 1),
            FeatureLocation(99498, 194508, 1),
            "tool", "T1PKS", 20, 20, "T1PKS",
        )
        created = creator([first, second, third], 4154245)

        assert len(created) == 2
        assert created[0].kind == created[0].kind.CHEMICAL_HYBRID
        assert set(created[0].protoclusters) == {first, second}
        assert created[1].kind == created[1].kind.SINGLE
        assert set(created[1].protoclusters) == {third}

    def test_merge(self):
        make = partial(create_cluster, 50, 60, 120, 170)
        # something that should not be merged
        isolated = create_cluster(0, 5, 10, 15, "isolated")
        # a pair of sets that should be merged, separated by the origin
        shared_ab = make("shared_ab")
        shared_ab.location = CompoundLocation([FeatureLocation(180, 200), FeatureLocation(0, 30)])
        a = create_cluster(15, 20, 40, 55, "a")  # pylint: disable=invalid-name
        b = create_cluster(170, 175, 185, 190, "b")  # pylint: disable=invalid-name
        # and an pair to merge that don't cross the origin
        c = make("c")  # pylint: disable=invalid-name
        cd = {c, make("d")}
        de = {c, make("e")}
        sets = [{a, shared_ab}, {isolated}, cd, de, {b, shared_ab}]
        results = merge(sets)
        results_as_sets = [set(clusters) for clusters in results]
        assert results_as_sets == [{a, b, shared_ab}, {isolated}, cd.union(de)]
