# License: GNU Affero General Public License v3 or later
# A copy of GNU AGPL v3 should have been included in this software package in LICENSE.txt.

""" The analysis sections of the TTA module.
    Checks for TTA codons within cluster features and adds a feature for each
    TTA codon found.

    Not recommended for low GC content sequences.
"""

import logging
from typing import Any, Dict

from antismash.common.secmet.feature import Feature, FeatureLocation
import antismash.common.module_results


class TTAResults(antismash.common.module_results.ModuleResults):
    """ Holds results for the TTA module by tracking locations of TTA codons."""
    schema_version = 1

    def __init__(self, record_id: str):
        super().__init__(record_id)
        self.codon_starts = []  # tuples of start and strand for each marker
        self.features = []  # features created for markers

    def new_feature_from_basics(self, start, strand) -> Feature:
        """ Constructs a new TTA marking feature from a start position and
            a strand
        """
        tta_feature = Feature(FeatureLocation(start, start + 3, strand), feature_type="misc_feature",
                              created_by_antismash=True)
        tta_feature.notes.append("tta leucine codon, possible target for bldA regulation")

        self.codon_starts.append((start, strand))
        self.features.append(tta_feature)

        return tta_feature

    def new_feature_from_other(self, feature, offset) -> Feature:
        """Create a misc_feature entry for a TTA codon on a given feature"""
        if feature.strand == 1:
            start = feature.location.start + offset
        else:
            start = feature.location.end - offset - 3

        return self.new_feature_from_basics(start, feature.strand)

    def to_json(self) -> Dict[str, Any]:
        """ Construct a JSON representation of this instance """
        starts = [{"start": start, "strand": strand} for start, strand in self.codon_starts]
        return {"TTA codons": starts,
                "schema_version": TTAResults.schema_version,
                "record_id": self.record_id}

    def add_to_record(self, record):
        """ Adds the found TTA features to the record """
        if record.id != self.record_id:
            raise ValueError("Record to store in and record analysed don't match")
        for feature in self.features:
            record.add_feature(feature)

    def __len__(self):
        return len(self.features)

    @staticmethod
    def from_json(json: Dict[str, Any], record) -> "TTAResults":
        """ Constructs a new TTAResults instance from a json format and the
            original record analysed.
        """
        if json["schema_version"] != TTAResults.schema_version:
            return None
        results = TTAResults(json["record_id"])
        for info in json["TTA codons"]:
            start = info["start"]
            strand = info["strand"]
            results.new_feature_from_basics(start, strand)
        return results


def detect(record, options) -> TTAResults:
    """Detect TTA codons"""
    assert options.tta
    logging.info("Detecting TTA codons")
    results = TTAResults(record.id)
    for feature in record.get_cds_features_within_clusters():
        sequence = feature.extract(record.seq)
        for i in range(0, len(sequence), 3):
            codon = sequence[i:i+3].lower()
            if codon == "tta":
                results.new_feature_from_other(feature, i)

    logging.debug("Detected %d TTA codons", len(results.features))

    return results
