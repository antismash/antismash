# License: GNU Affero General Public License v3 or later
# A copy of GNU AGPL v3 should have been included in this software package in LICENSE.txt.


"""
Reads a fasta file,
runs muscle over each sequence with KRdomains_muscle.fasta,
does some kind of Stachelhaus-like position fetch,
compares the generated signatures to a few string options to generate whether
        sequence is active and what stereochemistry it is
writes an outfile file in the tab-separated form: name, activity, stereochem
"""

from typing import Dict, Optional, Tuple

from antismash.common import path, subprocessing, utils
from antismash.modules.nrps_pks.data_structures import SimplePrediction, Prediction

_KR_DOMAINS_FILENAME = path.get_full_path(__file__, "data", "KRdomains_muscle.fasta")


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
    for name, seq in sorted(queries.items()):
        querysignames.append(name)
        muscle_dict = subprocessing.run_muscle_single(name, seq, _KR_DOMAINS_FILENAME)

        positions_act = [110, 134, 147, 151]  # active site
        positions_ste = [90, 91, 92, 139, 144, 147, 149, 151]  # stereochem

        refsequence = "MAPSI|PKS|CAM00062.1|Erythromycin_synthase_modules_1_and_2|Sacc_KR1"
        refseq = muscle_dict[refsequence]
        activity_signatures.append(utils.extract_by_reference_positions(muscle_dict[name], refseq, positions_act))
        stereochem_signatures.append(utils.extract_by_reference_positions(muscle_dict[name], refseq, positions_ste))

    # Check activity
    activity = {}  # type: Dict[str, Prediction]
    for name, signature in zip(querysignames, activity_signatures):
        if is_active(signature):
            activity[name] = SimplePrediction("kr_activity", "active")
        else:
            activity[name] = SimplePrediction("kr_activity", "inactive")

    # Predict stereochemistry
    stereochemistry = {}  # type: Dict[str, Prediction]
    for name, signature in zip(querysignames, stereochem_signatures):
        chem = predict_stereochemistry(signature)
        if chem:
            stereochemistry[name] = SimplePrediction("kr_stereochem", chem)

    return activity, stereochemistry
