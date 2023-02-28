# License: GNU Affero General Public License v3 or later
# A copy of GNU AGPL v3 should have been included in this software package in LICENSE.txt.

"""
Dynamic profile to predict cyanobactin precursors based on a shared leader peptide motif.

Based on doi: 10.1128/AEM.01061-09, the motif is M.KKN[IL].P....PV.R
Because in this manuscript the leading M is from the start codon, that would almost always
be conserved and thus doesn't carry a lot of information. We also don't restrict the pattern
to the start of a sequence, this allows detecting the microcyclamide precursors as well. If
if the motif starts later in these or if the start codons are misannotated
"""
from typing import Dict, List
import re


from antismash.common.hmm_rule_parser.structures import DynamicHit, DynamicProfile
from antismash.common.secmet import Record

# This is the name this profile will have in the cluster rules
PROFILE_NAME = "cyanobactin_precursor"
PROFILE_DESCRIPTION = "Pattern-based detection of cyanobactin precursor peptides"

ANCHOR = re.compile(r"(KK|PV)")  # Use both the KK and the PV motif to allow for mismatches
PATTERN: Dict[int, str] = {
    2: "K",
    3: "K",
    4: "N",
    5: "IL",
    7: "P",
    12: "P",
    13: "V",
    15: "R",
}

MAX_SCORE = len(PATTERN)  # Set dynamically so it updates if we ever need to change the pattern
ALLOWABLE_MISSES = 1
MIN_LEN = max(PATTERN) + 1  # If sequences aren't at least as long as the pattern, this won't work
MAX_LEN = 150
LAST_ALLOWABLE_START = 50  # We expect the pattern in the leader


def find_hits(record: Record, allowable_misses: int = ALLOWABLE_MISSES) -> Dict[str, List[DynamicHit]]:
    """Find all CDSes where the pattern is found"""
    hits: Dict[str, List[DynamicHit]] = {}

    for cds in record.get_cds_features():
        if not MIN_LEN <= len(cds.translation) <= MAX_LEN:
            continue

        bitscore = 0.0
        for hit in ANCHOR.finditer(cds.translation, endpos=LAST_ALLOWABLE_START):
            idx = hit.start()
            if hit.group(0) == "PV":
                idx -= 12
                if idx < 0:
                    continue
            if idx + 15 >= len(cds.translation):
                continue
            bitscore = check_match(cds.translation, idx, allowable_misses)
            if bitscore:
                hits[cds.get_name()] = [DynamicHit(cds.get_name(), PROFILE_NAME, bitscore=bitscore)]
                break

    return hits


def check_match(translation: str, start_position: int, allowable_misses: int) -> float:
    """Check if the pattern with at most allowable_misses mismatches is found in a translation"""
    misses = 0
    for idx, allowed in PATTERN.items():
        if translation[idx + start_position] not in allowed:
            misses += 1
            if misses > allowable_misses:
                return 0.0
    return (MAX_SCORE - misses) / MAX_SCORE


profile = DynamicProfile(
    PROFILE_NAME,
    PROFILE_DESCRIPTION,
    find_hits
)
