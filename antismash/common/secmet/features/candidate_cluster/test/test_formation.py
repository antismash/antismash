# License: GNU Affero General Public License v3 or later
# A copy of GNU AGPL v3 should have been included in this software package in LICENSE.txt.

# for test files, silence irrelevant and noisy pylint warnings
# pylint: disable=no-self-use,protected-access,missing-docstring

import unittest

from antismash.common.secmet import FeatureLocation
from antismash.common.secmet.features import CandidateCluster
from antismash.common.secmet.features.candidate_cluster import (
    create_candidates_from_protoclusters as creator,
)

from .test_candidate_cluster import create_cds, create_cluster


class TestCreation(unittest.TestCase):
    def test_creation_empty(self):
        assert creator([]) == []

    def test_creation_single(self):
        created = creator([create_cluster(3, 8, 71, 76, 'a')])
        print(created)
        assert len(created) == 1
        assert created[0].location == FeatureLocation(3, 76)
        assert created[0].kind == CandidateCluster.kinds.SINGLE

    def test_creation_neighbours(self):
        cluster = create_cluster(3, 8, 71, 76, 'a')
        extra_cluster = create_cluster(50, 100, 120, 170, 'b')
        created = creator([cluster, extra_cluster])
        print(created)
        assert len(created) == 3
        expected_location = FeatureLocation(cluster.location.start,
                                            extra_cluster.location.end)
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
        assert candidate_cluster.location == FeatureLocation(3, 170)

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
        assert candidate_cluster.location == FeatureLocation(3, 170)

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
        assert created[0].location == FeatureLocation(3, 270)
        assert created[0].kind == CandidateCluster.kinds.NEIGHBOURING
        assert created[0].protoclusters == (cluster, neighbour_cluster, hybrid_cluster, overlap_cluster)

        assert created[1].location == FeatureLocation(3, 180)
        assert created[1].kind == CandidateCluster.kinds.INTERLEAVED
        assert created[1].protoclusters == (cluster, hybrid_cluster, overlap_cluster)

        assert created[2].location == FeatureLocation(3, 170)
        assert created[2].kind == CandidateCluster.kinds.CHEMICAL_HYBRID
        assert created[2].protoclusters == (cluster, hybrid_cluster)

        assert created[3].location == FeatureLocation(50, 270)
        assert created[3].kind == CandidateCluster.kinds.SINGLE
        assert created[3].protoclusters == (neighbour_cluster,)

        assert created[4].location == FeatureLocation(450, 600)
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
            assert created[0].location == FeatureLocation(3, 250)
            assert created[0].kind == CandidateCluster.kinds.INTERLEAVED
            assert created[0].protoclusters == tuple(sorted([cluster, hybrid, contained, overlapping]))

            assert created[1].location == FeatureLocation(3, 176)
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
        assert created[0].location == FeatureLocation(30, 430)
        assert created[0].core_location == FeatureLocation(60, 410)
        assert created[0].kind == CandidateCluster.kinds.INTERLEAVED
        assert created[0].protoclusters == tuple(first_hybrid_clusters + second_hybrid_clusters + [single])

        assert created[1].location == FeatureLocation(30, 180)
        assert created[1].protoclusters == tuple(first_hybrid_clusters)
        assert created[2].location == FeatureLocation(90, 310)
        assert created[2].protoclusters == tuple(second_hybrid_clusters)
        for cand in created[1:3]:
            assert cand.kind == CandidateCluster.kinds.CHEMICAL_HYBRID

        assert created[3].location == FeatureLocation(1000, 1600)
        assert created[3].kind == CandidateCluster.kinds.CHEMICAL_HYBRID

    def test_interleaving_order(self):
        clusters = [create_cluster(1000, 1100, 1400, 1500, "a"),
                    create_cluster(1050, 2000, 3000, 4000, "b"), # sorts second due to neighbouring
                    create_cluster(1100, 1200, 1500, 1600, "c")]
        assert sorted(clusters) == clusters
        created = creator(clusters)
        assert len(created) == 3
        assert created[0].kind == CandidateCluster.kinds.NEIGHBOURING
        assert created[0].location == FeatureLocation(1000, 4000)
        assert created[1].kind == CandidateCluster.kinds.INTERLEAVED
        assert created[1].location == FeatureLocation(1000, 1600)
        assert created[2].kind == CandidateCluster.kinds.SINGLE
        assert created[2].location == FeatureLocation(1050, 4000)

    def test_contained_neighbours(self):
        big = create_cluster(0, 100, 130, 230, "a")
        small = create_cluster(10, 20, 30, 40, "b")

        created = creator([big, small])

        assert len(created) == 2
        assert created[0].kind == CandidateCluster.kinds.NEIGHBOURING
        assert created[0].location == big.location

        assert created[1].kind == CandidateCluster.kinds.SINGLE
        assert created[1].location == small.location

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
