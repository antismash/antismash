# License: GNU Affero General Public License v3 or later
# A copy of GNU AGPL v3 should have been included in this software package in LICENSE.txt.

"""
Dynamic profile to predict darobactin radical SAM enzymes in the vicinity of the precursor.

Based on doi: 10.1038/s41586-019-1791-1, the motif is W.W.K..
"""

from antismash.common.hmm_rule_parser.structures import (
    DynamicHit,
    DynamicProfile,
    ProfileHit,
)
from antismash.common.secmet import Record
from antismash.detection.hmm_detection.dynamic_profiles.darobactin_precursor import ANCHOR as MOTIF, MIN_LEN, MAX_LEN
from antismash.detection.hmm_detection.dynamic_profiles._utils import find_motif_around_anchor

# This is the name this profile will have in the cluster rules
PROFILE_NAME = "darobactin_rSAM"
PROFILE_DESCRIPTION = "Detection of darobactin-related rSAM based on precursor presence"

MAX_DIST = 10000

def find_hits(record: Record, hmmer_hits: dict[str, list[ProfileHit]]) -> dict[str, list[DynamicHit]]:
    """Find all CDSes where the pattern is found"""
    hits: dict[str, list[DynamicHit]] = {}

    for anchor in record.get_cds_features():
        found_rsam = False
        found_spasm = False

        for hit in hmmer_hits.get(anchor.get_name(), []):
            if hit.query_id == "PF04055":
                found_rsam = True
            elif hit.query_id in ("TIGR04085", "SPASM"):
                found_spasm = True
        if not all((found_spasm, found_rsam)):
            continue

        precursors = find_motif_around_anchor(
            record=record,
            anchor=anchor,
            motif=MOTIF,
            max_dist=MAX_DIST,
            min_len=MIN_LEN,
            max_len=MAX_LEN,
            early_abort=True,
        )

        if precursors:
            hits[anchor.get_name()] = [DynamicHit(anchor.get_name(), PROFILE_NAME, bitscore=1.0)]

    return hits


profile = DynamicProfile(
    PROFILE_NAME,
    PROFILE_DESCRIPTION,
    find_hits
)
