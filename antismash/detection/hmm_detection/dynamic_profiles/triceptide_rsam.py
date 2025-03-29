# License: GNU Affero General Public License v3 or later
# A copy of GNU AGPL v3 should have been included in this software package in LICENSE.txt.

"""
Dynamic profile to predict triceptide radical SAM enzymes in the vicinity of the precursor.

Based on doi: 10.1021/jacs.2c00521, the motif is AxYxDxP, HAASL, or YxRxHxxHxR
"""
import re

from antismash.common.all_orfs import find_all_orfs, get_trimmed_orf
from antismash.common.hmm_rule_parser.structures import (
    DynamicHit,
    DynamicProfile,
    ProfileHit,
)
from antismash.common.secmet import Record
from antismash.common.secmet.features import CDSCollection
from antismash.common.secmet.locations import FeatureLocation

# This is the name this profile will have in the cluster rules
PROFILE_NAME = "triceptide_rSAM"
PROFILE_DESCRIPTION = "Detection of triceptide-related rSAM based on precursor presence"

MOTIF = re.compile(r"(A.Y.D.P|HAASL|Y.R.H..H.R|WDN)")
MAX_DIST = 1000

MIN_LEN = 10

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
                break
            elif hit.query_id in ("PF04055"):
                found_rsam = True
            elif hit.query_id in ("TIGR04085", "SPASM"):
                found_spasm = True
        if not all((found_spasm, found_rsam)):
            continue
        begin = max(0, int(anchor.location.start) - MAX_DIST)
        end = min(len(record), int(anchor.location.end) + MAX_DIST)

        area = CDSCollection(location=FeatureLocation(begin, end), feature_type="Area")
        existing_features = record.get_cds_features_within_location(area.location)
        found = False
        for cds in existing_features:
            if MOTIF.search(cds.translation) is not None:
                hits[anchor.get_name()] = [DynamicHit(anchor.get_name(), PROFILE_NAME, bitscore=1.0)]
                found = True
                break

        # No need to search for a precursor
        if found:
            continue

        new_features = find_all_orfs(record, area, min_length=MIN_LEN)
        for cds in new_features:
            new_orf = get_trimmed_orf(cds, record, min_length=MIN_LEN, label=cds.get_name())
            if new_orf:
                cds = new_orf
            if MOTIF.search(cds.translation) is not None:
                hits[anchor.get_name()] = [DynamicHit(anchor.get_name(), PROFILE_NAME, bitscore=1.0)]
                found = True
                break

    return hits


profile = DynamicProfile(
    PROFILE_NAME,
    PROFILE_DESCRIPTION,
    find_hits
)
