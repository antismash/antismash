# License: GNU Affero General Public License v3 or later
# A copy of GNU AGPL v3 should have been included in this software package in LICENSE.txt.

""" Predicts CAL domain specificities by Minowa et al. method"""

from typing import Dict

from antismash.common import path
from .base import run_minowa, Prediction

CAL_DOMAINS_PATH = path.get_full_path(__file__, "data", "CAL_domains_muscle.fasta")


def run_minowa_cal(sequence_info: Dict[str, str]) -> Dict[str, Prediction]:
    """ Predicts CAL domain specificities by Minowa et al. method"""
    return run_minowa(sequence_info=sequence_info,
                      startpos=43,
                      muscle_ref=CAL_DOMAINS_PATH,
                      ref_sequence="Q54297_CAL1",
                      positions_file=path.get_full_path(__file__, "data", "CALpositions.txt"),
                      data_dir=path.get_full_path(__file__, "data", 'CAL_HMMs'),
                      hmm_names=["Acetyl-CoA",
                                 "AHBA",
                                 "fatty_acid",
                                 "NH2",
                                 "shikimic_acid"])
