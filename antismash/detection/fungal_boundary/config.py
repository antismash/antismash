# License: GNU Affero General Public License v3 or later
# A copy of GNU AGPL v3 should have been included in this software
# package in LICENSE.txt.

"""Parameters and derived constants for fungal boundary prediction."""

UPSTREAM = 1000                 # promoter bases before the gene start
DOWNSTREAM = 50                 # promoter bases after the gene start
LINEAR_UPSTREAM = 300           # linear motifs only

ARM = 3                         # bases per repeat arm
MAX_SPACER = 15                 # max bases allowed between the two arms
LINEAR_LENGTH = 6               # bases in a linear motif
MAX_COUNTS = 3                  # max hits of one motif counted per promoter

N_FLANKING = 120
MIN_GENES = 2                   # smallest window tested
MAX_GENES = 50                  # largest window tested
MIN_GENES_WITH_MOTIFS = 3       # min promoters in a window with motif

NUM_TEST_ANCHORS = 200          # random genes used for calibration
SEED = 1                        # keeps the random genes reproducible
MIN_FLANK = 20                  # genes required on each side of a random gene

P_VAL_MAX = 0.05                # empirical p-value a prediction must beat
MAX_PREDICTED = 30              # largest window still accepted
MIN_LLR_PER_GENE = 0.75         # score a window must reach per gene
MIN_LINEAR_SUPPORT = 0.75       # fraction of window with a linear motif
BUFFER_GENES = 2                # genes added in the accepted window

MODES = ("anchor", "adjacent", "all")
TYPES = ("inverted", "direct", "linear")

SPACERS = MAX_SPACER + 1
# motif ID slots to reserve for each repeat type

# 4 ** ARM: possible DNA sequences for one arm. With four bases
REPEATS = 4 ** ARM * SPACERS
DIRECT_START = REPEATS
LINEAR_START = 2 * REPEATS
N_MOTIFS = LINEAR_START + 4 ** LINEAR_LENGTH
MOTIF_TYPE_RANGE = {
    "inverted": (0, DIRECT_START),
    "direct": (DIRECT_START, LINEAR_START),
    "linear": (LINEAR_START, N_MOTIFS),
}
