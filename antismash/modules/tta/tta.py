# License: GNU Affero General Public License v3 or later
# A copy of GNU AGPL v3 should have been included in this software package in LICENSE.txt.

""" The analysis sections of the TTA module.
    Checks for TTA codons within cluster features and adds a feature for each
    TTA codon found.

    Not recommended for low GC content sequences.
"""

import logging
from typing import Any, Dict, List, Tuple, Optional

from antismash.common.secmet import Record
from antismash.common.secmet.features import Feature, FeatureLocation
from antismash.common.module_results import ModuleResults
from antismash.config import ConfigType, get_config

Codon = Tuple[int, int]


class TTAResults(ModuleResults):
    """ Holds results for the TTA module by tracking locations of TTA codons."""
    schema_version = 2

    def __init__(self, record_id: str, gc_content: float, threshold: float) -> None:
        super().__init__(record_id)
        self.codon_starts: List[Codon] = []  # tuples of start and strand for each marker
        self.features: List[Feature] = []  # features created for markers
        self.gc_content = float(gc_content)
        self.threshold = float(threshold)

    def new_feature_from_basics(self, start: int, strand: int) -> Feature:
        """ Constructs a new TTA marking feature from a start position and
            a strand
        """
        tta_feature = Feature(FeatureLocation(start, start + 3, strand), feature_type="misc_feature",
                              created_by_antismash=True)
        tta_feature.notes.append("tta leucine codon, possible target for bldA regulation")

        self.codon_starts.append((start, strand))
        self.features.append(tta_feature)

        return tta_feature

    def new_feature_from_other(self, feature: Feature, offset: int) -> Feature:
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
                "record_id": self.record_id,
                "gc_content": self.gc_content,
                "threshold": self.threshold}

    def add_to_record(self, record: Record) -> None:
        """ Adds the found TTA features to the record """
        if record.id != self.record_id:
            raise ValueError("Record to store in and record analysed don't match")
        for feature in self.features:
            record.add_feature(feature)

    def __len__(self) -> int:
        return len(self.features)

    @staticmethod
    def from_json(json: Dict[str, Any], record: Record) -> Optional["TTAResults"]:
        """ Constructs a new TTAResults instance from a json format and the
            original record analysed.
        """
        if json["schema_version"] != TTAResults.schema_version:
            return None

        options = get_config()
        results = TTAResults(json["record_id"], json["gc_content"], options.tta_threshold)
        # if old results were excluding based on too low a GC content, rerun
        if json["threshold"] > results.gc_content and options.tta_threshold <= results.gc_content:
            return None
        # otherwise, if the threshold is now too high, skip all the codons
        if json["gc_content"] >= get_config().tta_threshold:
            for info in json["TTA codons"]:
                start = info["start"]
                strand = info["strand"]
                results.new_feature_from_basics(start, strand)
        return results


def detect(record: Record, options: ConfigType) -> TTAResults:
    """ Find TTA codons in a record

        Arguments:
            record: the record to search
            options: an antismash config object

        Returns:
            a TTAResults instance with all detected TTA codons
    """
    gc_content = record.get_gc_content()
    results = TTAResults(record.id, gc_content, options.tta_threshold)

    if gc_content < options.tta_threshold:
        logging.info("Skipping TTA codon detection, GC content too low: %2.0f%%", gc_content * 100)
        return results

    logging.info("Detecting TTA codons")
    logging.debug("GC content: %2.0f%%", gc_content * 100)

    for feature in record.get_cds_features_within_regions():
        sequence = feature.extract(record.seq)
        for i in range(0, len(sequence), 3):
            codon = sequence[i:i+3].lower()
            if codon == "tta":
                results.new_feature_from_other(feature, i)

    logging.debug("Detected %d TTA codons", len(results.features))

    return results
