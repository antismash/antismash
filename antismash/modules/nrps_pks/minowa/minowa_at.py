# License: GNU Affero General Public License v3 or later
# A copy of GNU AGPL v3 should have been included in this software package in LICENSE.txt.

""" Predicts AT domain specificities by Minowa et al. method """

from typing import Dict

from antismash.common import path
from antismash.modules.nrps_pks.minowa.base import run_minowa, Prediction

AT_DOMAINS_PATH = path.get_full_path(__file__, "data", "AT_domains_muscle.fasta")


def run_minowa_at(sequence_info: Dict[str, str]) -> Dict[str, Prediction]:
    """ Predicts AT domain specificities by Minowa et al. method """
    return run_minowa(sequence_info=sequence_info,
                      startpos=7,
                      muscle_ref=AT_DOMAINS_PATH,
                      ref_sequence="P0AAI9_AT1",
                      positions_file=path.get_full_path(__file__, "data", "ATpositions.txt"),
                      data_dir=path.get_full_path(__file__, "data", 'AT_HMMs'),
                      hmm_names=["2-Methylbutyryl-CoA",
                                 "Acetyl-CoA",
                                 "CHC-CoA",
                                 "fatty_acid",
                                 "Isobutyryl-CoA",
                                 "Methoxymalonyl-CoA",
                                 "Propionyl-CoA",
                                 "3-Methylbutyryl-CoA",
                                 "Benzoyl-CoA",
                                 "Ethylmalonyl-CoA",
                                 "inactive",
                                 "Malonyl-CoA",
                                 "Methylmalonyl-CoA",
                                 "trans-1,2-CPDA"])
