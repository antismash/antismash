# License: GNU Affero General Public License v3 or later
# A copy of GNU AGPL v3 should have been included in this software package in LICENSE.txt.

"""
Dynamic profile to predict atropopeptide precursors based on the conserved recognition sequence.

Based on doi: 10.1039/d4sc03469d, the motif is KPLK, with some variations
"""
import re


from antismash.common.hmm_rule_parser.structures import (
    DynamicHit,
    DynamicProfile,
    ProfileHit,
)
from antismash.common.secmet import Record

# This is the name this profile will have in the cluster rules
PROFILE_NAME = "atropopeptide_precursor"
PROFILE_DESCRIPTION = "Pattern-based detection of atropopeptide precursor peptides"

ANCHOR = re.compile(r"(?:[KREQP]S|KP)LK")

MIN_LEN = 12  # The final compound core is 5-6 AAs, the recognition motif ends 2 AAs before the core
MAX_LEN = 80
MAX_DIST = 3000


def find_hits(record: Record, hmmer_hits: dict[str, list[ProfileHit]]) -> dict[str, list[DynamicHit]]:
    """Find all CDSes where the pattern is found"""
    hits: dict[str, list[DynamicHit]] = {}

    for p450 in record.get_cds_features():
        found = False
        for hit in hmmer_hits.get(p450.get_name(), []):
            if hit.query_id in ("p450", "atropopeptide_p450"):
                found = True
                break
        if not found:
            continue

        location = record.extend_location(p450.location, MAX_DIST)
        existing_features = record.get_cds_features_within_location(location)
        for cds in existing_features:
            if not MIN_LEN <= len(cds.translation) <= MAX_LEN:
                continue

            if ANCHOR.search(cds.translation) is not None:
                hits[cds.get_name()] = [DynamicHit(cds.get_name(), PROFILE_NAME, bitscore=1.0)]

    return hits


profile = DynamicProfile(
    PROFILE_NAME,
    PROFILE_DESCRIPTION,
    find_hits
)
