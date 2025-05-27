# License: GNU Affero General Public License v3 or later
# A copy of GNU AGPL v3 should have been included in this software package in LICENSE.txt.

"""
Functions used by multiple dynamic profiles

"""
import re

from antismash.common.all_orfs import find_all_orfs, get_trimmed_orf
from antismash.common.secmet import Record
from antismash.common.secmet.features import CDSCollection
from antismash.common.secmet.features.cds_feature import CDSFeature

def find_motif_around_anchor(
        record: Record,
        anchor: CDSFeature,
        motif: re.Pattern,
        *,
        max_dist: int,
        min_len: int,
        max_len: int = 0,
        early_abort: bool = False,
    ) -> list[CDSFeature]:
    """ Find a sequence motif in small ORFs around the anchor and return the found ORFs

        Arguments:
            record: the Record the hits are located in
            anchor: the CDSFeature around which to look for small ORFs with the motif
            motif: the compiled regex Pattern to look for in the ORFs' translations
            max_dist: how many kb to extend around the anchor
            min_len: the minimum size of the small ORF in AAs
            max_len: the maximum size of the small ORF in AAs. if 0, then don't impose a limit
            early_abort: when True, abort after the first CDS was found, otherwise, find all the ORFs

        Returns:
            a list of ORFs around the anchor that contain the motif

    """
    found: list[CDSFeature] = []

    if max_len < 0:
        raise ValueError("max_len needs to be >= 0")

    search_location = record.extend_location(anchor.location, max_dist)
    area = CDSCollection(location=search_location, feature_type="Area")
    existing_features = record.get_cds_features_within_location(area.location)
    for cds in existing_features:
        if len(cds.translation) < min_len:
            continue
        if max_len and len(cds.translation) > max_len:
            continue
        if motif.search(cds.translation) is not None:
            found.append(cds)
            if early_abort:
                break

    # No need to continue with gene finding if we found something and should abort early
    if found and early_abort:
        return found

    # min_len is in AAs
    min_nucleotides = min_len * 3
    new_features = find_all_orfs(record, area, min_length=min_nucleotides)
    for cds in new_features:
        new_orf = get_trimmed_orf(cds, record, min_length=min_nucleotides, label=cds.get_name())
        if new_orf:
            cds = new_orf
        if max_len and len(cds.translation) > max_len:
            continue
        if motif.search(cds.translation) is not None:
            found.append(cds)
            if early_abort:
                break

    return found
