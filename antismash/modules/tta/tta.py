# License: GNU Affero General Public License v3 or later
# A copy of GNU AGPL v3 should have been included in this software package in LICENSE.txt.

import logging

from antismash.common import deprecated

def detect(seq_record, options):
    """Detect TTA codons"""
    assert options.tta
    logging.info("Detecting TTA codons")
     #TODO: change to new secmet structures
    cds_features = deprecated.get_withincluster_cds_features(seq_record)
    new_features = []
    for feature in cds_features:
        sequence = feature.extract(seq_record.seq)
        for i in range(0, len(sequence), 3):
            codon = sequence[i:i+3].lower()
            if codon == "tta":
                tta = _create_tta_feature(feature, i)
                seq_record.features.append(tta)
                new_features.append(tta)
    logging.debug("Detected %d TTA codons", len(new_features))
    return new_features


def _create_tta_feature(feature, offset):
    """Create a misc_feature entry for a TTA codon on a given feature"""
    if feature.strand == 1:
        start = feature.location.start + offset
        end = start + 3
    else:
        end = feature.location.end - offset
        start = end - 3

    loc = deprecated.FeatureLocation(start, end, feature.strand)

    qualifiers = {
        "note": ["tta leucine codon, possible target for bldA regulation"],
        "tool": ["antiSMASH"]
    }
    tta_feature = deprecated.SeqFeature(loc, type="misc_feature", qualifiers=qualifiers)
    return tta_feature
