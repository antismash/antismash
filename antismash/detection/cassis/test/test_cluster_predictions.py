# License: GNU Affero General Public License v3 or later
# A copy of GNU AGPL v3 should have been included in this software package in LICENSE.txt.

# for test files, silence irrelevant and noisy pylint warnings
# pylint: disable=use-implicit-booleaness-not-comparison,protected-access,missing-docstring

import unittest

from antismash.common.secmet.features import Gene, FeatureLocation
from antismash.detection.cassis.cluster_prediction import check_cluster_predictions, \
                    ClusterMarker, ClusterPrediction, create_predictions
from antismash.detection.cassis.islands import Island
from antismash.detection.cassis.motifs import Motif
from antismash.detection.cassis.promoters import Promoter, CombinedPromoter

from .test_cassis import create_fake_record


class TestClusterPredictions(unittest.TestCase):
    def test_check_cluster_predictions(self):
        seq_record = create_fake_record()
        promoters = [Promoter("gene1", 1, 5),
                     Promoter("gene2", 6, 10),
                     CombinedPromoter("gene3", "gene4", 11, 15)]
        ignored_genes = [  # see captured logging
            Gene(FeatureLocation(1, 5), locus_tag="gene5")
        ]
        clusters = [ClusterPrediction(ClusterMarker("gene1", Motif(3, 3, score=1)),
                                      ClusterMarker("gene4", Motif(3, 3, score=1)))]
        expected = [ClusterPrediction(ClusterMarker("gene1", Motif(3, 3, score=1)),
                                      ClusterMarker("gene4", Motif(3, 3, score=1)))]
        expected[0].start.promoter = "gene1"
        expected[0].end.promoter = "gene3+gene4"
        expected[0].genes = 4
        expected[0].promoters = 3

        assert check_cluster_predictions(clusters, seq_record, promoters, ignored_genes) == expected
        # -------------------- >> begin captured logging << --------------------
        # root: INFO: Best prediction (most abundant): 'gene1' -- 'gene4'
        # root: INFO: Upstream cluster border located at or next to sequence
        #               record border, prediction could have been truncated by record border
        # root: INFO: Downstream cluster border located at or next to sequence
        #               record border, prediction could have been truncated by record border
        # root: INFO: Gene 'gene2' is part of the predicted cluster, but it is
        #               overlapping with another gene and was ignored
        # root: INFO: Gene 'gene2' could have effected the cluster prediction
        # --------------------- >> end captured logging << ---------------------

    def test_sort_by_abundance(self):
        islands = []

        # island 1: [gene1 -- gene2]
        motif = Motif(0, 3, score=3, hits={"gene1": 1, "gene2": 1})
        islands.append(Island(Promoter("gene1", 1, 1), Promoter("gene2", 2, 2), motif))
        # island 2: [gene2 -- gene5]
        motif = Motif(3, 0, score=2, hits={"gene2": 1, "gene3": 1, "gene4": 1, "gene5": 1})
        islands.append(Island(Promoter("gene2", 2, 2), Promoter("gene5", 5, 5), motif))
        # island 3: [gene1 -- gene5]
        motif = Motif(3, 3, score=1, hits={"gene1": 1, "gene2": 1, "gene3": 1, "gene4": 1, "gene5": 1})
        islands.append(Island(Promoter("gene1", 1, 1), Promoter("gene5", 5, 5), motif))

        # left border: 2x gene1, 1x gene2
        # right border: 2x gene5, 1x gene2

        expected_clusters = []
        # cluster 1: [gene1 -- gene5] --> abundance 2+2 (most abundant)
        start = ClusterMarker("gene1", Motif(3, 3, score=1))
        start.abundance = 2
        end = ClusterMarker("gene5", Motif(3, 3, score=1))
        end.abundance = 2
        expected_clusters.append(ClusterPrediction(start, end))
        # cluster 3: [gene2 -- gene5] --> abundance 1+2, score 2+1 (better/lower)
        start = ClusterMarker("gene2", Motif(3, 0, score=2))
        start.abundance = 1
        end = ClusterMarker("gene5", Motif(3, 3, score=1))
        end.abundance = 2
        expected_clusters.append(ClusterPrediction(start, end))
        # cluster 2: [gene1 -- gene2] --> abundance 2+1, score 1+3 (worse, higher)
        start = ClusterMarker("gene1", Motif(3, 3, score=1))
        start.abundance = 2
        end = ClusterMarker("gene2", Motif(0, 3, score=3))
        end.abundance = 1
        expected_clusters.append(ClusterPrediction(start, end))
        # cluster 4: [gene2 -- gene2] --> abundance 1+1
        start = ClusterMarker("gene2", Motif(3, 0, score=2))
        start.abundance = 1
        end = ClusterMarker("gene2", Motif(0, 3, score=3))
        end.abundance = 1
        expected_clusters.append(ClusterPrediction(start, end))
        # abundance: as high as possible
        # score: as low as possible

        self.assertEqual(create_predictions(islands), expected_clusters)
