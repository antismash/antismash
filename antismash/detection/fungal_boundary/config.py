# License: GNU Affero General Public License v3 or later
# A copy of GNU AGPL v3 should have been included in this software
# package in LICENSE.txt.

"""Parameters and motif ID ranges for fungal boundary analysis."""

MIN_GENES_PER_CONTIG = 1000      # minimum CDS features for boundary analysis
MIN_GENES_WITH_MOTIF = 3         # minimum genes containing a motif for enrichment
MAX_GENES_IN_BGC = 30           # maximum genes in a scored window
UPSTREAM = 1000                 # promoter bases before the gene start
DOWNSTREAM = 50                 # promoter bases after the gene start
LINEAR_UPSTREAM = 300           # linear motifs only

ARM = 3                         # bases per repeat arm
MAX_SPACER = 15                 # max bases allowed between the two arms
LINEAR_LENGTH = 6               # bases in a linear motif
MAX_COUNTS = 3                  # max hits of one motif counted per promoter
MAX_RANDOM_ANCHORS = 100
MIN_RANDOM_ANCHORS = 10
MAX_RANDOM_ANCHOR_OVERLAP = 28
RANDOM_SEED = 1
MIN_LLR = 20                    # minimum LLR required for empirical p-value testing
P_VAL_MAX = 0.05                # empirical p-value a prediction must beat

SPACERS = MAX_SPACER + 1

# Number of motif IDs reserved separately for inverted and direct repeats.
MOTIF_IDS_PER_REPEAT_TYPE = 4 ** ARM * SPACERS
DIRECT_MOTIF_ID_START = MOTIF_IDS_PER_REPEAT_TYPE
LINEAR_MOTIF_ID_START = 2 * MOTIF_IDS_PER_REPEAT_TYPE
N_MOTIFS = LINEAR_MOTIF_ID_START + 4 ** LINEAR_LENGTH

# A motif_id's numeric range identifies its type:
#   inverted repeat: 0 <= motif_id < DIRECT_MOTIF_ID_START
#   direct repeat:   DIRECT_MOTIF_ID_START <= motif_id < LINEAR_MOTIF_ID_START
#   linear motif:    LINEAR_MOTIF_ID_START <= motif_id < N_MOTIFS
# With searchsorted(..., side="right"), these boundaries give type indices
# 0 (inverted), 1 (direct), and 2 (linear).
MOTIF_TYPE_ID_BOUNDARIES = (DIRECT_MOTIF_ID_START, LINEAR_MOTIF_ID_START)
