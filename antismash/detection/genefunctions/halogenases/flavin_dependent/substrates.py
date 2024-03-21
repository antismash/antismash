import os
import pathlib
from typing import List, Union, Optional

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
GENERAL_FDH_PROFILES = [HmmSignature("all_general_FDH",
                                     "Member of the Flavin-dependent halogenase family",
                            100, get_full_path(pathlib.Path(__file__).parents[1], "data","halogenases",
                                               "all_general_FDH.hmm")),
                        HmmSignature("unconventional_FDH",
                                     "Unconventional member of the Flavin-dependent halogenase family",
                            100, get_full_path(pathlib.Path(__file__).parents[1], "data", "halogenases",
                                               "unconventional_FDH.hmm"))]

GENERAL_FDH_MOTIFS = {"W.W.I.": [206, 207, 208, 209, 210, 211],
                        "F.*P.*S.G": [280, 281, 282, 283, 284, 285, 286,
                                        287, 288, 289, 290, 291, 292, 293,
                                        294, 295, 296, 297, 298, 299, 300,
                                        301, 302, 303, 304]
                                        }
