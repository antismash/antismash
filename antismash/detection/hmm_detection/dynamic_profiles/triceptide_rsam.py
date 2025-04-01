# License: GNU Affero General Public License v3 or later
# A copy of GNU AGPL v3 should have been included in this software package in LICENSE.txt.

"""
Dynamic profile to predict triceptide radical SAM enzymes in the vicinity of the precursor.

Based on doi: 10.1021/jacs.2c00521, the motif is AxYxDxP, HAASL, or YxRxHxxHxR
Based on doi: 10.1021/acschembio.2c00621, WDN is another motif
"""

from antismash.common.hmm_rule_parser.structures import (
    DynamicHit,
    DynamicProfile,
    ProfileHit,
)
from antismash.common.secmet import Record
from antismash.detection.hmm_detection.dynamic_profiles.triceptide_precursor import ANCHOR as MOTIF, MIN_LEN
from antismash.detection.hmm_detection.dynamic_profiles._utils import find_motif_around_anchor

# This is the name this profile will have in the cluster rules
PROFILE_NAME = "triceptide_rSAM"
PROFILE_DESCRIPTION = "Detection of triceptide-related rSAM based on precursor presence"

MAX_DIST = 1000


def find_hits(record: Record, hmmer_hits: dict[str, list[ProfileHit]]) -> dict[str, list[DynamicHit]]:
    """Find all CDSes where the pattern is found"""
    hits: dict[str, list[DynamicHit]] = {}

    for anchor in record.get_cds_features():
        found_rsam = False
        found_spasm = False

        for hit in hmmer_hits.get(anchor.get_name(), []):
            if hit.query_id in ("SAM_SPASM_FxsB", "rSAM_GlyRichRpt", "rSAM_XyeB", "PTHR43273"):
                found_rsam = True
                found_spasm = True
            found_rsam |= hit.query_id in ("PF04055")
            found_spasm |= hit.query_id in ("TIGR04085", "SPASM")
            if found_rsam and found_spasm:
                break
        if not all((found_spasm, found_rsam)):
            continue

        precursors = find_motif_around_anchor(
            record=record,
            anchor=anchor,
            motif=MOTIF,
            max_dist=MAX_DIST,
            min_len=MIN_LEN,
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
