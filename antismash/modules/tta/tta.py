# License: GNU Affero General Public License v3 or later
# A copy of GNU AGPL v3 should have been included in this software package in LICENSE.txt.

import json
import logging

from antismash.common.secmet.feature import Feature
from antismash.common import deprecated
import antismash.common.module_results

class TTAResults(antismash.common.module_results.ModuleResults):
    schema_version = 1
    def __init__(self, record_id):
        super().__init__(record_id)
        self.codon_starts = []
        self.features = []

    def new_feature_from_basics(self, start, strand):
        loc = deprecated.FeatureLocation(start, start + 3, strand)

        qualifiers = {
            "note": ["tta leucine codon, possible target for bldA regulation"],
            "tool": ["antiSMASH"]
        }
        #TODO change to using a Feature directly
        tta_feature = Feature.from_biopython(deprecated.SeqFeature(loc, type="misc_feature", qualifiers=qualifiers))

        self.codon_starts.append((start, strand))
        self.features.append(tta_feature)

        return tta_feature

    def new_feature_from_other(self, feature, offset):
        """Create a misc_feature entry for a TTA codon on a given feature"""
        if feature.strand == 1:
            start = feature.location.start + offset
        else:
            start = feature.location.end - offset - 3

        return self.new_feature_from_basics(start, feature.strand)

    def to_json(self):
        starts = [{"start" : start, "strand" : strand} for start, strand in self.codon_starts]
        return {"TTA codons" : starts,
                "schema_version" : TTAResults.schema_version,
                "record_id" : self.record_id}

    def add_to_record(self, record):
        if record.id != self.record_id:
            raise ValueError("Record to store in and record analysed don't match")
        for feature in self.features:
            record.add_feature(feature)

    def __len__(self):
        return len(self.features)

    @staticmethod
    def from_json(data):
        if isinstance(data, str):
            data = json.loads(data)
        if data["schema_version"] != TTAResults.schema_version:
            return None
        results = TTAResults(data["record_id"])
        for info in data["TTA codons"]:
            start = info["start"]
            strand = info["strand"]
            results.new_feature_from_basics(start, strand)
        return results

def create_results_from_json(data):
    if data["schema_version"] != TTAResults.schema_version: # or <, whatever is possible
        raise ValueError("Result schema version mismatch, cannot parse")
    results = TTAResults(data["record_id"])

    for codon in data["TTA codons"]:
        start = codon["start"]
        strand = codon["strand"]
        results.new_feature_from_basics(start, strand)
    return results

def detect(seq_record, options):
    """Detect TTA codons"""
    assert options.tta
    logging.info("Detecting TTA codons")
    results = TTAResults(seq_record.id)
    for feature in seq_record.get_cds_features_within_clusters():
        sequence = feature.extract(seq_record.seq)
        for i in range(0, len(sequence), 3):
            codon = sequence[i:i+3].lower()
            if codon == "tta":
                results.new_feature_from_other(feature, i)

    logging.debug("Detected %d TTA codons", len(results.features))

    results.add_to_record(seq_record)

    return results
