# License: GNU Affero General Public License v3 or later
# A copy of GNU AGPL v3 should have been included in this software package in LICENSE.txt.

""" A very basic genefinding system. Finds all ORFs in a sequence, with some
    extra conditions applied.

    Not intended for genefinding an unannotated sequence, but for finding extra
    ORFs that have been skipped in an annotated sequence such as RiPP precursors.
"""


from typing import Iterable, List, Optional, Tuple

from Bio.SeqFeature import FeatureLocation

from antismash.common.secmet import CDSFeature, Record
from antismash.common.secmet.features import CDSCollection

START_CODONS = ('ATG', 'GTG', 'TTG')
STOP_CODONS = ('TAA', 'TAG', 'TGA')


def scan_orfs(seq: str, direction: int, offset: int = 0, minimum_length: int = 60
              ) -> List[FeatureLocation]:
    """ Scan for open reading frames on a given sequence.
        Skips all ORFs with a size less than the given minimum nucleotides.

        Arguments:
            seq: the sequence to examine
            direction: the search direction to use (all ORFs will use this as the strand)
            offset: an offset to add to any location discovered
            minimum_length: the minimum length, in nucleotides, for ORFs

        Returns:
            a list of FeatureLocations for each ORF, ordered by ascending position
    """
    seq = seq.upper()

    matches = []
    # cache the sequence length
    seq_len = len(seq)
    for frame in [0, 1, 2]:
        start = None
        for i in range(frame, seq_len - 2, 3):
            codon = seq[i:i+3]
            # use the earliest possible start
            if start is None and codon in START_CODONS:
                start = i
                continue
            if codon in STOP_CODONS:
                # skip stops without a matching start
                if start is None:
                    continue
                end = i + 2  # include the full codon
                # cull genes smaller than the cutoff
                if end - start < minimum_length:
                    start = None
                    continue
                # finally, calculate the appropriate location on the record
                if direction == 1:
                    loc_start = start + offset
                    loc_end = end + offset + 1
                else:
                    loc_start = seq_len + offset - end - 1
                    loc_end = seq_len + offset - start
                matches.append(FeatureLocation(loc_start, loc_end, direction))
                start = None
    return sorted(matches, key=lambda x: min(x.start, x.end))


def create_feature_from_location(record: Record, location: FeatureLocation,
                                 label: Optional[str] = None) -> CDSFeature:
    """ Creates a CDS feature covering the provided location.

        Arguments:
            record: The Record the CDSFeature will belong to, used to generate
                    the feature translation
            location: The FeatureLocation specifying the location of the CDSFeature
            label: The locus tag, protein id, and gene name to use for the new
                   CDSFeature

        Returns:
            The CDSFeature created.
    """
    if label is None:
        digits = len(str(len(record)))
        label = 'allorf_{start:0{digits}}_{end:0{digits}}'.format(
            digits=digits, start=(location.start + 1), end=location.end
        )
    feature = CDSFeature(location, str(record.get_aa_translation_from_location(location)),
                         locus_tag=label, protein_id=label, gene=label)
    feature.created_by_antismash = True
    return feature


def find_intergenic_areas(start: int, end: int, cds_features: Iterable[CDSFeature],
                          min_length: int = 0, padding: int = 0) -> List[Tuple[int, int]]:
    """ Finds intergenic areas between the given CDS features.

        Arguments:
            start: the minimum coordinate to use
            end: the maximum coordinate to use
            cds_features: the collection of CDS features
                            (assumed to be ordered by at least start location)
            min_length: the minimum length of intergenic area to keep
            padding: the allowable overlap with existing CDS features

        Returns:
            a list of tuples, each indicating the start and end of an intergenic area
    """
    intergenic_areas = []
    last = start
    for cds in cds_features:
        # if there's a gap of sufficient size, add it
        if cds.location.start > last:
            intergenic_areas.append((max(start, last - padding),
                                     min(end, int(cds.location.start) + padding)))
            last = int(cds.location.end)
            continue
        # in case of existing CDS features overlapping, update the last end position
        if cds.location.start <= last <= cds.location.end:
            last = int(cds.location.end)
    if last < end:
        intergenic_areas.append((max(start, last - padding), end))
    return list(filter(lambda area: area[1] - area[0] >= min_length, intergenic_areas))


def find_all_orfs(record: Record, area: Optional[CDSCollection] = None,
                  min_length: int = 60, max_overlap: int = 10) -> List[CDSFeature]:
    """ Find ORFs within intergenic areas of the given record or subset of the record.

        Can (and should) be limited to just within a specific section of the record.

        Arguments:
            record: the record to search
            area: the specific CDSCollection to search within, or None
            min_length: the minimum length of ORFs to report, in nucleotides
            max_overlap: the maximum allowable bases of overlap with existing CDS features

        Returns:
            a list of CDSFeatures, one for each ORF
    """
    # Get sequence for the range
    offset = 0
    seq = record.seq
    existing: Iterable[CDSFeature] = record.get_cds_features()
    if area:
        seq = area.extract(seq)
        offset = area.location.start
        existing = record.get_cds_features_within_location(area.location,
                                                           with_overlapping=True)
    intergenic_areas = find_intergenic_areas(offset, offset + len(seq), existing,
                                             min_length=min_length, padding=max_overlap)

    # Find orfs throughout the range
    locations = []
    for start, end in intergenic_areas:
        chunk = seq[start:end]
        locations.extend(scan_orfs(chunk, 1, start))
        locations.extend(scan_orfs(chunk.reverse_complement(), -1, start))

    new_features = []
    for location in locations:
        new_features.append(create_feature_from_location(record, location))

    return sorted(new_features)
