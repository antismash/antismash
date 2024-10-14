# License: GNU Affero General Public License v3 or later
# A copy of GNU AGPL v3 should have been included in this software package in LICENSE.txt.

"""
Dynamic profile to predict triceptide precursors based on the cross-linked amino acids.

Based on doi: 10.1021/jacs.2c00521, the motif is AxYxDxP, HAASL, or YxRxHxxHxR
"""
from typing import Dict, List
import re


from antismash.common.hmm_rule_parser.structures import DynamicHit, DynamicProfile
from antismash.common.secmet import Record

# This is the name this profile will have in the cluster rules
PROFILE_NAME = "triceptide_precursor"
PROFILE_DESCRIPTION = "Pattern-based detection of triceptide precursor peptides"

ANCHOR = re.compile(r"(A.Y.D.P|HAASL|Y.R.H..H.R)")

MIN_LEN = 10
MAX_LEN = 100

def find_hits(record: Record) -> Dict[str, List[DynamicHit]]:
    """Find all CDSes where the pattern is found"""
    hits: Dict[str, List[DynamicHit]] = {}

    for cds in record.get_cds_features():
        if not MIN_LEN <= len(cds.translation) <= MAX_LEN:
            continue

        if ANCHOR.search(cds.translation) is not None:
            hits[cds.get_name()] = [DynamicHit(cds.get_name(), PROFILE_NAME, bitscore=10.0)]

    return hits


profile = DynamicProfile(
    PROFILE_NAME,
    PROFILE_DESCRIPTION,
    find_hits
)
