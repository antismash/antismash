# License: GNU Affero General Public License v3 or later
# A copy of GNU AGPL v3 should have been included in this software package in LICENSE.txt.

""" Classes and helpers for candidate cluster features

CandidateClusters contain one or more Protocluster features.

There are four kinds of candidate cluster:
    - chemical hybrid:
        contains clusters which share Cluster-defining CDSFeatures, will also
        include clusters within that shared range that do not share a CDS provided
        that they are completely contained within the candidate cluster border,
        e.g.
            ---##A###############C#---   <- Cluster 1 with definition CDSes A and C
             --##A##--                   <- Cluster 2 with definition CDS A
                      --#B#--            <- Cluster 3 with definition CDS B
                              --#C#--    <- Cluster 4 with definition CDS C
                                   -#D#- <- Cluster 5 with definition CDS D
            Since clusters 1 and 2 share a CDS that defines those clusters, a
            chemical hybrid candidate cluster exists. Clusters 1 and 4 also share a
            defining CDS, so the hybrid candidate cluster now contains clusters 1, 2 and 4.

            Cluster 3 does not share a defining CDS with either cluster 1, 2 or 4,
            but because it is interleaved into a chemical hybrid it is included
            under the assumption that it is relevant to the other clusters.
      NOTE: This may change so that there is also an 'interleaved' candidate cluster,
            where the only difference between the interleaved and the hybrid is
            that Cluster 3 would contribute it's product to the interleaved and
            not the hybrid.
    - interleaved:
        contains clusters which do not share Cluster-defining CDS features, but
        their core locations overlap,
        e.g.
            ---#A###A###A---      <- Cluster 1 with defining CDSes marked A
               ---B##B####B---    <- Cluster 2 with defining CDSes marked B
                      ---C###C--- <- Cluster 3 with defining CDSes marked C
            Since none of the clusters share any defining CDS with any other cluster,
            it is not a chemical hybrid. All three clusters would be part of an
            interleaved candidate cluster, since A overlaps with B and B overlaps with C.
    - neighbouring:
        contains clusters which transitively overlap in their neighbourhoods
        (the '-' sections in the examples above). In the chemical hybrid example,
        as all clusters overlap in some way, all 5 would be part of a neighbouring
        candidate cluster (with clusters 1-4 also being part of a hybrid candidate cluster).
        Every cluster in a 'neighbouring' cluster will also belong to one of the
        other kinds of candidate cluster.
    - single:
        the kind for all candidate clusters where only one cluster is contained,
        only exists for consistency of access. A 'single' candidate cluster will not
        exist for a cluster which is contained in either a chemical hybrid or
        an interleaved candidate cluster. In the chemical hybrid example, only cluster 5
        would be in a 'single' candidate cluster as well as in the 'neighbouring' candidate cluster

"""

from .formation import create_candidates_from_protoclusters
from .structures import CandidateCluster, CandidateClusterKind
