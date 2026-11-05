# License: GNU Affero General Public License v3 or later
# A copy of GNU AGPL v3 should have been included in this software package in LICENSE.txt.

"""
Module to detect transcription factor binding sites (TFBSs)
using position weight matrices (PWMs) in BGCs.
"""

from dataclasses import dataclass
from enum import IntEnum, auto
import itertools
import logging

from typing import Any, Dict, Iterator, List, Optional, Tuple

from collections import Counter
from Bio.Seq import Seq
from MOODS import tools, scan

from antismash.config import get_config
from antismash.common import json as jsonlib, path
from antismash.common.module_results import ModuleResults
from antismash.common.secmet import Record
from antismash.common.secmet.features import Feature, Region
from antismash.common.secmet.locations import (
    CompoundLocation,
    FeatureLocation,
    Location,
)


PWM_PATH = path.get_full_path(__file__, 'data', 'PWMs.json')


class Confidence(IntEnum):
    """ Defined and sortable values for confidence terms """
    WEAK = auto()
    MEDIUM = auto()
    STRONG = auto()

    def __str__(self) -> str:
        return self.name.lower()


@dataclass
class Matrix:
    """ Class to store position weight matrices (Matrix) """
    name: str
    pwm: List[List[float]]
    max_score: float
    min_score: float
    description: str
    species: str
    link: str
    consensus: str
    _threshold: float = -1.

    @property
    def score_threshold(self) -> float:
        """ The minimum threshold required for a hit's score to be considered good """
        if self._threshold < 0:
            self._threshold = (self.min_score + self.max_score) / 2
        return self._threshold

    def get_score_confidence(self, score: float) -> Confidence:
        """ Returns the confidence level for the given score """
        if score <= self.min_score:
            return Confidence.WEAK
        if score >= self.score_threshold:
            return Confidence.STRONG
        return Confidence.MEDIUM

    def to_json(self) -> Dict[str, Any]:
        """ Returns a JSON-ready representation of the instance """
        return {key: val for key, val in vars(self).items() if not key.startswith("_")}

    @staticmethod
    def from_json(name: str, data: Dict[str, Any]) -> "Matrix":
        """ Reconstructs an instance from a JSON representation """
        return Matrix(name, **data)


@dataclass
class Hit:
    """ Class to store hits retrieved from the analysis tool MOODS """
    pos: int
    score: float
    strand: int

    @classmethod
    def from_moods_match(cls, match: scan.match, strand: int, offset: int) -> "Hit":
        """ Builds an instance from a MOODS result """
        if strand not in [-1, 1]:
            raise ValueError(f"invalid strand: {strand}")
        return Hit(match.pos + offset, match.score, strand)


@dataclass
class TFBSHit:
    """ Class to store the transcription factor binding site hits """
    name: str
    start: int
    species: str
    link: str
    description: str
    consensus: str
    confidence: Confidence
    strand: int
    score: float
    max_score: float

    def to_json(self) -> Dict[str, Any]:
        """ Returns a JSON-ready representation of the instance """
        data = dict(vars(self))
        data["confidence"] = str(data["confidence"])
        return data

    @staticmethod
    def from_json(data: Dict[str, Any]) -> "TFBSHit":
        """ Reconstructs an instance from a JSON representation """
        assert isinstance(data["confidence"], str)
        data = dict(data)
        data["confidence"] = Confidence[data["confidence"].upper()]
        return TFBSHit(**data)

    @classmethod
    def from_matrix_hit(cls, hit: Hit, matrix: Matrix) -> "TFBSHit":
        """ Constructs a detailed hit instance from some hit details and matrix information

            Arguments:
                hit: the hit
                matrix: the matrix used to find the hit

            Returns:
                a new instance of TFBSHit
        """
        if hit.strand not in [-1, 1]:
            raise ValueError(f"invalid strand: {hit.strand}")
        return cls(matrix.name, hit.pos, matrix.species, matrix.link, matrix.description, matrix.consensus,
                   matrix.get_score_confidence(hit.score), hit.strand, hit.score, matrix.max_score)

    def __lt__(self, other: Any) -> bool:
        if not isinstance(other, type(self)):
            raise TypeError(f"cannot order {type(self)} and {type(other)}")
        if self.confidence > other.confidence:  # by default, sort by decreasing confidence
            return True
        if self.start < other.start:  # increasing start
            return True
        return self.name < other.name  # and alphabetically

    def __len__(self) -> int:
        return len(self.consensus)


class TFBSFinderResults(ModuleResults):
    """ Results class for the TFBS finder analysis """
    schema_version = 1

    def __init__(self, record_id: str, pvalue: float, start_overlap: int,
                 upstream_range: int, hits_by_region: Dict[int, List[TFBSHit]]) -> None:
        super().__init__(record_id)
        self.pvalue = pvalue  # the p-value used for threshold setting within MOODS
        self.start_overlap = start_overlap
        self.upstream_range = upstream_range
        self.hits_by_region = hits_by_region

    def get_hits_by_region(self, region_number: int, confidence: Confidence = None,
                           allow_better: bool = False) -> List[TFBSHit]:
        """ Returns hits from the given region, optionally limiting by confidence

            Arguments:
                region_number: the region number to pull results for
                confidence: a confidence to limits results to
                allow_better: if true, also includes all confidences stronger than the given confidence

            Returns:
                a list of the filtered hits
        """
        if allow_better and not confidence:
            raise ValueError("cannot specify 'allow_better' without also specifying 'confidence'")
        hits = self.hits_by_region.get(region_number, [])
        if allow_better and confidence:
            # < since decreasing order
            return sorted([hit for hit in hits if hit.confidence >= confidence])
        if confidence:
            return sorted([hit for hit in hits if hit.confidence == confidence])
        return sorted(hits)

    def to_json(self) -> Dict[str, Any]:
        """ Constructs a JSON representation of this instance """
        hits_by_region = {}
        for region, hits in self.hits_by_region.items():
            hits_by_region[region] = [hit.to_json() for hit in hits]

        return {
            "schema_version": self.schema_version,
            "pvalue": self.pvalue,
            "start_overlap": self.start_overlap,
            "upstream_range": self.upstream_range,
            "record_id": self.record_id,
            "hits_by_region": hits_by_region,
        }

    def add_to_record(self, record: Record) -> None:
        """ Adds the analysis results to the record """
        if record.id != self.record_id:
            raise ValueError("Record to store in and record analysed don't match")
        for hits in self.hits_by_region.values():
            for hit in hits:
                start = hit.start % len(record)
                end = (hit.start + len(hit.consensus)) % len(record)
                if start > end:
                    location = CompoundLocation([
                        FeatureLocation(start, len(record), hit.strand),
                        FeatureLocation(0, end, hit.strand),
                    ])
                else:
                    location = FeatureLocation(start, end, hit.strand)
                tfbs_feature = Feature(location,
                                       feature_type="misc_feature", created_by_antismash=True)
                tfbs_feature.notes.append(f"TFBS match to {hit.name}, {hit.description}, "
                                          f"confidence: {hit.confidence}, "
                                          f"score: {round(hit.score, 2)}")
                record.add_feature(tfbs_feature)

    @staticmethod
    def from_json(json: Dict[str, Any], record: Record) -> Optional["TFBSFinderResults"]:
        """ Constructs a new results instance from a JSON format and the
            original record analysed.
        """
        # check that the previous data version is the same as current, if not, discard the results
        if json["schema_version"] != TFBSFinderResults.schema_version:
            return None

        if record.id != json.get("record_id"):
            logging.warning("TFBS results are for different record, discarding previous results")
            return None

        prev_pval = json["pvalue"]
        assert isinstance(prev_pval, float)

        prev_range = json["start_overlap"]
        assert isinstance(prev_range, int)

        prev_upstream = json["upstream_range"]
        assert isinstance(prev_upstream, int)

        # fetch the cutoffs
        options = get_config()
        upstream_range, start_overlap = parse_range(options.tfbs_range)
        if options.tfbs_pvalue != prev_pval:
            logging.debug("TFBSFinderResults p-value has changed, discarding previous results")
            return None

        if start_overlap != prev_range or upstream_range != prev_upstream:
            logging.debug("TFBSFinderResults search range has changed, discarding previous results")
            return None

        hits_by_region = {}

        for name, hits in json["hits_by_region"].items():
            hits_by_region[int(name)] = [TFBSHit.from_json(hit) for hit in hits]

        return TFBSFinderResults(json["record_id"], options.tfbs_pvalue,
                                 start_overlap, upstream_range, hits_by_region)


def parse_range(value: str) -> Tuple[int, int]:
    """ Converts a range option into the distances it refers to, upstream and
        downstream of a gene start position

        Arguments:
            value: the option value, in the form "upstream:downstream"

        Returns:
            a tuple of the upstream and downstream distances
    """
    upstream, separator, downstream = value.partition(":")
    if not separator:
        raise ValueError(f"expected an 'upstream:downstream' range, not: {value}")
    return int(upstream), int(downstream)


def get_sequence_gc_content(sequence: str) -> float:
    """ Calculates the proportion of the given nucleotide sequence that are G or C

        Arguments:
            sequence: the sequence to calculate

        Returns:
            a float in the range 0 to 1
    """
    counter = Counter(sequence)
    return sum(counter[char] for char in "GC") / len(sequence)


def get_bg_distribution(sequence: Seq) -> Tuple[float, float, float, float]:
    """ Calculates the background distribution of a sequence based on GC percentage
        in a format that MOODS will accept

        Arguments:
            sequence: the sequence to calculate the background distribution of

        Returns:
            a tuple of AT-, GT-, GT-, and AT-proportions of the sequence
    """
    gc_percentage = get_sequence_gc_content(sequence) / 2
    at_percentage = 0.5 - gc_percentage
    # this repetition is intentional, as that's what MOODS needs
    return at_percentage, gc_percentage, gc_percentage, at_percentage


def load_matrices(pwm_file: str) -> List[Matrix]:
    """ Loads the position weight matrices from the matrix file

        Arguments:
            pwm_file: the path to the matrix file

        Returns:
            a list of Matrix instances, one for each matrix in the file
    """
    matrices = []
    with open(pwm_file, "r", encoding="utf-8") as file:
        regulators = jsonlib.load(file)
        for reg_id, reg_info in regulators.items():
            reg_info["min_score"] = reg_info.get("min_score", 0)
            matrices.append(Matrix(reg_id, **reg_info))
    return matrices


def run_moods(sequence: Seq, background: Tuple[float, float, float, float],
              pwms: List[Matrix], pvalue: float, offset: int) -> List[List[Hit]]:
    """ Runs MOODS over the given sequence, with the given matrices

        Arguments:
            sequence: the sequence to scan
            background: the background distribution matrix
            pwms: the matrices to scan with
            pvalue: the minimum threshold for MOODS hits
            offset: the offset to add to hit positions

        Returns:
            a list of lists of hits, one sublist for each matrix
    """
    matrices = [pwm_info.pwm for pwm_info in pwms]
    thresholds = [tools.threshold_from_p(m, background, pvalue) for m in matrices]
    forward = scan.scan_dna(str(sequence), matrices, background, thresholds)
    reverse = scan.scan_dna(str(sequence.reverse_complement()), matrices, background, thresholds)
    results = []
    for i, matrix in enumerate(matrices):
        starts: dict[int, Hit] = {}
        for raw_hit in forward[i]:
            hit = Hit.from_moods_match(raw_hit, 1, offset)
            existing = starts.get(hit.pos)
            if existing and existing.score > hit.score:
                continue
            starts[hit.pos] = hit
        for raw_hit in reverse[i]:
            start = len(sequence) - raw_hit.pos - len(matrix[0]) + offset
            hit = Hit(start, raw_hit.score, -1)
            existing = starts.get(hit.pos)
            if existing and existing.score > hit.score:
                continue
            starts[hit.pos] = hit
        results.append(list(starts.values()))
    return results


def are_co_transcribed(upstream: Location, downstream: Location) -> bool:
    """ Identify likely co-transcribed genes based on their intergenic distances.

        Arguments:
            upstream: the location with the lower coordinates
            downstream: the location with the higher coordinates

        Returns:
            True if the two locations are predicted to be part of the same transcript
    """
    if upstream.strand != downstream.strand:
        return False
    distance = downstream.start - upstream.end
    return -90 <= distance <= 170


def get_valid_areas(start: int, end: int, locations: Iterator[Location],
                    start_overlap: int, upstream_range: int) -> list[tuple[int, int]]:
    """ Finds areas within the section where a binding site is considered valid.
        Generally these are intergenic areas, but some overlaps with genes are possible.

        Arguments:
            start: the start coordinate of the section to search in
            end: the end coordinate of the section to search in
            locations: the existing locations to excise
            start_overlap: the allowable overlap when the start of two locations overlap
                           and are on reverse strands

        Returns:
            a list of tuples, each being a start and end coordinate of an area,
                ordered by increasing start coordinate
    """
    original_end = end
    areas: List[Tuple[int, int]] = []

    def add_area(start: int, end: int) -> None:
        if start < end:
            # making sure to strip BeforePosition and other variants from locations
            areas.append((int(start), int(end)))

    # add initial range
    try:
        initial = next(locations)
    except StopIteration:
        # no locations to excise, so return the whole block
        return [(start, end)]
    if initial.strand == 1:
        add_area(max(start, initial.start - upstream_range), initial.start + start_overlap)

    previous = initial
    for current in locations:
        if current.start < previous.end:  # overlaps with a previous CDS
            # overlaps and same strand or of two ends can be ignored
            # but overlaps of two CDS starts are interesting
            if current.strand == 1 and previous.strand == -1:
                # with no gap between them, both regions are added in full
                add_area(current.start - upstream_range, current.start + start_overlap)
                add_area(previous.end - start_overlap, previous.end + upstream_range)

            if current.end > previous.end:
                previous = current
            continue

        # otherwise, skip gaps between the ends of two genes
        if current.strand == -1 and previous.strand == 1:
            previous = current
            continue

        # also skip when genes are predicted to be co-transcribed
        if are_co_transcribed(previous, current):
            previous = current
            continue

        if previous.strand == -1:
            add_area(previous.end - start_overlap, min(previous.end + upstream_range, current.start))

        if current.strand == 1:
            add_area(max(current.start - upstream_range, previous.end), current.start + start_overlap)

        previous = current

    # add final range, but only if the last gene's start codon faces it
    if previous.strand == -1:
        add_area(previous.end - start_overlap, min(original_end, previous.end + upstream_range))

    # the regions of two genes can overlap, so merge any that do
    areas.sort()
    merged: List[Tuple[int, int]] = []
    for area in areas:
        if merged and area[0] <= merged[-1][1]:
            merged[-1] = (merged[-1][0], max(merged[-1][1], area[1]))
            continue
        merged.append(area)
    return merged


def filter_hits(matrices: List[Matrix], areas: List[Tuple[int, int]],
                moods_results: List[List[Hit]]) -> List[TFBSHit]:
    """ Filters hits to those in valid areas within the region.

        Arguments:
            matrices: a list of matrices from which the hits were generated
            areas: valid areas that hits can be within
            moods_results: the results from moods, each sublist matching the same

        Returns:
            a list of hits that fall in valid areas
    """
    if len(matrices) != len(moods_results):
        raise ValueError("")
    if not areas:
        return []
    results = []
    for matrix, hits in zip(matrices, moods_results):
        if not hits:
            continue
        hits = sorted(hits, key=lambda hit: hit.pos)
        hit_index = 0
        area_index = 0
        while hit_index < len(hits) and area_index < len(areas):
            hit = hits[hit_index]
            area = areas[area_index]

            if hit.pos < area[0]:  # skip hits before the current area
                hit_index += 1
                continue

            if area[1] <= hit.pos:  # skip areas before the current hit
                area_index += 1
                continue

            assert area[0] <= hit.pos < area[1], f"{hit.pos}, {area}"
            hit_index += 1
            if matrix.get_score_confidence(hit.score) < Confidence.MEDIUM:  # exclude weak matches
                continue
            results.append(TFBSHit.from_matrix_hit(hit, matrix))
    return results


def get_cross_origin_areas(record: Record, region: Region, start_overlap: int,
                           upstream_range: int) -> list[tuple[int, int]]:
    """ Finds areas within the section where a binding site is considered valid.
        Generally these are intergenic areas, but some overlaps with genes are possible.

        Arguments:
            record: the record containing the region
            region: the region to find areas in
            start_overlap: the allowable overlap with a start codon
            upstream_range: the distance before a start codon that a binding site can be found in

        Returns:
            a list of tuples, each being a start and end coordinate of an area,
                ordered by increasing start coordinate
    """
    pre_origin = region.cds_children.pre_origin
    cross_origin = region.cds_children.cross_origin
    post_origin = region.cds_children.post_origin

    record_length = len(record)
    locations: list[Location] = []

    for cds in pre_origin:
        locations.append(FeatureLocation(int(cds.start), int(cds.end), cds.location.strand))

    for cds in cross_origin:
        parts = cds.location.parts
        # the parts are ordered by transcription, so reverse strand genes start at the end
        pre_part, post_part = (parts[-1], parts[0]) if cds.location.strand == -1 else (parts[0], parts[-1])
        locations.append(FeatureLocation(int(pre_part.start), int(post_part.end) + record_length,
                                         cds.location.strand))

    for cds in post_origin:
        locations.append(FeatureLocation(int(cds.start) + record_length, int(cds.end) + record_length,
                                         cds.location.strand))

    # to match the linear version, ensure all locations are in order of lowest start coordinate
    locations.sort(key=lambda location: (location.start, location.end))

    return get_valid_areas(int(region.start), int(region.end) + record_length, iter(locations),
                           start_overlap, upstream_range)


def run_tfbs_finder(record: Record, pvalue: float, start_overlap: int, upstream_range: int,
                    matrix_path: str = PWM_PATH) -> TFBSFinderResults:
    """Run TFBS finder on a given record

        Arguments:
            record: the record to search
            pvalue: the cutoff p-value
            start_overlap: the upstream search range
            matrix_path: the path to a file containing position weight matrices

        Returns:
            a TFBSResults instance with all detected binding sites
    """
    hits_by_region = {}
    matrices = load_matrices(matrix_path)
    for region in record.get_regions():
        sequence = region.extract(record.seq)
        background = get_bg_distribution(sequence)
        moods_results = run_moods(sequence, background, matrices, pvalue, region.start)
        areas: list[tuple[int, int]] = []
        if region.crosses_origin():
            areas = get_cross_origin_areas(record, region, start_overlap, upstream_range)
        else:
            locations = (cds.location for cds in region.cds_children)
            areas = get_valid_areas(region.start, region.end, locations, start_overlap, upstream_range)
        hits_by_region[region.get_region_number()] = filter_hits(matrices, areas, moods_results)
    return TFBSFinderResults(record.id, pvalue, start_overlap, upstream_range, hits_by_region)
