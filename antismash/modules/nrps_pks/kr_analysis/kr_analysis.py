# License: GNU Affero General Public License v3 or later
# A copy of GNU AGPL v3 should have been included in this software package in LICENSE.txt.

""" Analysis of PKS KR domain stereochemistry and activity
"""

from typing import Dict, Optional, Tuple

from antismash.common import path, utils
from antismash.common.brawn import get_cached_alignment, get_aligned_pair
from antismash.modules.nrps_pks.data_structures import SimplePrediction, Prediction

DATA_DIR = path.get_full_path(__file__, "data")
KR_DOMAINS_PATH = path.get_full_path(__file__, "data", "KRdomains_muscle.fasta")


def is_active(signature: str) -> bool:
    """ Returns True if a KR domain is active based on its signature """
    if signature[0] == "E" and signature[1] in "SAG" and signature[2] == "H" and signature[3] == "H":
        return True
    if signature[0] == "K" and signature[1] in "SAG" and signature[2] == "Y" and signature[3] in "NG":
        return True
    return False


def predict_stereochemistry(signature: str) -> Optional[str]:
    """ Predicts stereochemistry of a KR domain from its signature """
    stereochemistry = None
    if signature[0:3] != "LDD" and signature[3] == "W" and signature[5:] == "YAN":
        if signature[4] != "H":
            stereochemistry = "A1"
        else:
            stereochemistry = "A2"
    elif signature[0:3] == "LDD" and signature[5] == "Y" and signature[7] == "N":
        if signature[6] != "P":
            stereochemistry = "B1"
        else:
            stereochemistry = "B2"
    elif signature[5] != "Y":
        stereochemistry = "C1"
    elif signature[7] != "N":
        stereochemistry = "C2"
    return stereochemistry


def run_kr_analysis(queries: Dict[str, str]) -> Tuple[Dict[str, Prediction], Dict[str, Prediction]]:
    """ Extract activity and stereochemistry signatures from KR domains

        Arguments:
            queries: a mapping of query CDS name to sequence

        Returns:
            a pair of dicts, one mapping query name to activity bool,
                the other mapping query name to stereochemistry (e.g. A2)
    """
    querysignames = []
    activity_signatures = []
    stereochem_signatures = []
    alignment = get_cached_alignment(KR_DOMAINS_PATH, DATA_DIR)
    reference = "MAPSI|PKS|CAM00062.1|Erythromycin_synthase_modules_1_and_2|Sacc_KR1"
    for name, seq in sorted(queries.items()):
        querysignames.append(name)
        aligned, refseq = get_aligned_pair(seq, reference, alignment)
        positions_act = [110, 134, 147, 151]  # active site
        positions_ste = [90, 91, 92, 139, 144, 147, 149, 151]  # stereochem

        activity_sig = utils.extract_by_reference_positions(aligned, refseq, positions_act)
        assert activity_sig
        activity_signatures.append(activity_sig)
        stereochem_sig = utils.extract_by_reference_positions(aligned, refseq, positions_ste)
        assert stereochem_sig
        stereochem_signatures.append(stereochem_sig)

    # Check activity
    activity: Dict[str, Prediction] = {}
    for name, signature in zip(querysignames, activity_signatures):
        if is_active(signature):
            activity[name] = SimplePrediction("kr_activity", "active")
        else:
            activity[name] = SimplePrediction("kr_activity", "inactive")

    # Predict stereochemistry
    stereochemistry: Dict[str, Prediction] = {}
    for name, signature in zip(querysignames, stereochem_signatures):
        chem = predict_stereochemistry(signature)
        if chem:
            stereochemistry[name] = SimplePrediction("kr_stereochem", chem)

    return activity, stereochemistry
