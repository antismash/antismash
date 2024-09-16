# License: GNU Affero General Public License v3 or later
# A copy of GNU AGPL v3 should have been included in this software package in LICENSE.txt.

"""
Dynamic profile to predict atropopeptide p450s on top of the generic p450 pHMM.

This is done by searching the vicinity of the p450 for the presence of a small ORF
with an atropopeptide precursor recognition sequence.

Based on doi: 10.1002/anie.202208361
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
PROFILE_NAME = "atropopeptide_p450"
PROFILE_DESCRIPTION = "Precursor recognition sequence-based detection of atropopeptide p450 enzymes"


MAX_DIST = 3000
MOTIF = re.compile(r"(?:[KREQP]S|KP)LK")
MIN_LEN = 40  # The final compound core is 5-6 AAs, and the recognition motif ends 2 AAs before the core



def find_hits(record: Record, hmmer_hits: dict[str, list[ProfileHit]]) -> dict[str, list[DynamicHit]]:
    """Find all CDSe that are p450s and have a small ORF with the pattern close by"""
    hits: dict[str, list[DynamicHit]] = {}

    for anchor in record.get_cds_features():
        # if it's not a p450, skip
        found = False
        for hit in hmmer_hits.get(anchor.get_name(), []):
            if hit.query_id == "p450":
                found = True
                break
        if not found:
            continue

        begin = max(0, int(anchor.location.start) - MAX_DIST)
        end = min(len(record), int(anchor.location.end) + MAX_DIST)

        area = CDSCollection(location=FeatureLocation(begin, end), feature_type="Area")
        existing_features = record.get_cds_features_within_location(area.location)
        found = False
        for cds in existing_features:
        # find existing features using record.get_cds_features_within_location(area.location)

            if MOTIF.search(cds.translation) is not None:
                hits[anchor.get_name()] = [DynamicHit(anchor.get_name(), PROFILE_NAME, bitscore=1.0)]
                found = True
                break

        # No need to continue until we actually add the new pecursors to the record.
        if found:
            continue

        # find unannotated features using all_orfs.find_all_orfs(record, area, min_size=MIN_SIZE)
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

