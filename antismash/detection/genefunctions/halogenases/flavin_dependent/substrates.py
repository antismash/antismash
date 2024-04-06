# License: GNU Affero General Public License v3 or later
# A copy of GNU AGPL v3 should have been included in this software package in LICENSE.txt.

# for test files, silence irrelevant and noisy pylint warnings
# pylint: disable=use-implicit-booleaness-not-comparison,protected-access,missing-docstring

from pathlib import Path

from antismash.common.path import get_full_path
from antismash.common.signature import HmmSignature

SHORT_DESCRIPTION = """Categorization of halogenases based on family and function.
                       It utilizes pHMMs with thresholds and the presence or absence of
                       characteristic/signature residues to identify the substrate or
                       the number of halogenation.
                       Signature residues for the pyrrole substrate were determined by
                       using information from an experimental study (https://doi.org/10.1073/pnas.1519695113).
                       Signature residues for all the other substrates were extracted computationally,
                       by looking at the positions that have the same value.
                    """

GENERAL_FDH_PROFILES = [HmmSignature("all_conventional_FDH",
                                     "Member of the Flavin-dependent halogenase family",
                                     100, get_full_path(str(Path(__file__).parents[0]),
                                                        "data", "all_conventional_FDH.hmm")),
                        HmmSignature("unconventional_FDH",
                                     "Unconventional flavin-dependent halogenase",
                                     100, get_full_path(str(Path(__file__).parents[0]),
                                                        "data", "unconventional_FDH.hmm"))]

ALL_FDH_PROFILES = get_full_path(str(Path(__file__).parents[0]),
                                 "data", "FDH.hmm")

GENERAL_FDH_MOTIFS = {"W.W.I.": [206, 207, 208, 209, 210, 211],
                      "F.*P.*S.G": [280, 281, 282, 283, 284, 285, 286,
                                    287, 288, 289, 290, 291, 292, 293,
                                    294, 295, 296, 297, 298, 299, 300,
                                    301, 302, 303, 304]
                                        }
