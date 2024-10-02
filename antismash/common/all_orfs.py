# License: GNU Affero General Public License v3 or later
# A copy of GNU AGPL v3 should have been included in this software package in LICENSE.txt.

""" A very basic genefinding system. Finds all ORFs in a sequence, with some
    extra conditions applied.

    Not intended for genefinding an unannotated sequence, but for finding extra
    ORFs that have been skipped in an annotated sequence such as RiPP precursors.
"""

import re
from typing import Iterable, List, Optional, Tuple

from Bio.SeqFeature import FeatureLocation

from antismash.common.secmet import CDSFeature, Record
from antismash.common.secmet.features import CDSCollection

START_CODONS = ('ATG', 'GTG', 'TTG')
STOP_CODONS = ('TAA', 'TAG', 'TGA')
RIBOSOMAL_BINDING_SITE = re.compile('(?:GGA[AG]G)')
FLEXIBLE_RIBOSOMAL_BINDING_SITE = re.compile('(?=(?:[^GA]*[GA]){6})[GACT]{8}')


def get_trimmed_orf(orf: CDSFeature, record: Record, include: int = None,
                    min_length: int = 0, max_length: int = None,
                    label: Optional[str] = None) -> Optional[CDSFeature]:
    """ Creates a trimmed ORF to start at the latest possible start codon in the
        given ORF that satisfies the given options

        Arguments:
            orf: the CDS feature to create a trimmed ORF from
            record: the parent Record of the CDS
            include: a coordinate within the ORF's DNA sequence to include
            min_length: the minimum length to trim to
            max_length: the maximum length to trim to,
                        if not given, then it defaults to the length of the given ORF
            label: the locus tag for the newly created CDSFeature, if given

        Returns:
            a new CDSFeature with a smaller location or None if no shorter
            version could be found
    """
    seq = orf.extract(record.seq)
    # set defaults
    if max_length is None:
        max_length = len(seq)
    if include is None:
        include = len(seq)

    # check for bad values
    if min_length > max_length:
        raise ValueError("minimum length cannot be greater than maximum length")
    if min_length > len(seq):  # min length and sequence length are incompatible
        return None
    if max_length < len(seq) - include:  # max length and include are incompatible
        return None

    # construct the search range, while ensuring that codons are in the same
    # frame as the original
    start = max(0, len(seq) - (max_length - (max_length % 3)))
    end = min(len(seq) - min_length, include)

    # gather possible alternative start coordinates
    starts = []
    for i in range(start, end, 3):
        if seq[i:i + 3] not in START_CODONS:
            continue
        assert min_length <= len(seq) - i <= max_length
        starts.append(i)
    # if no shorter version can be found, don't create a new feature
    if not starts:
        return None

    # otherwise, use the start which gives the smallest possible ORF
    if orf.location.strand == 1:
        start = orf.location.start + starts[-1]
        end = orf.location.end
    else:
        start = orf.location.start
        end = orf.location.end - starts[-1]

    location = FeatureLocation(start, end, orf.location.strand)
    return create_feature_from_location(record, location, label=label)


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
                                 label: Optional[str] = None, rbs: Optional[bool] = None) -> CDSFeature:
    """ Creates a CDS feature covering the provided location.

        Arguments:
            record: The Record the CDSFeature will belong to, used to generate
                    the feature translation
            location: The FeatureLocation specifying the location of the CDSFeature
            label: The locus tag, protein id, and gene name to use for the new
                   CDSFeature
            rbs: The boolean if a ribosomal binding site can found upstream of the CDS

        Returns:
            The CDSFeature created.
    """
    if label is None:
        digits = len(str(len(record)))
        label = 'allorf_{start:0{digits}}_{end:0{digits}}'.format(
            digits=digits, start=(location.start + 1), end=location.end
        )
    translation = str(record.get_aa_translation_from_location(location))
    # always start with methionine for CDS features, even if it had an alternate start codon
    if translation[0] != "M":
        translation = "M" + translation[1:]
    feature = CDSFeature(location, translation,
                         locus_tag=label, protein_id=label, gene=label, rbs=rbs)
    feature.created_by_antismash = True
    return feature


def has_rbs(sequence: str, flexible: bool = False) -> bool:
    """ Determines if ribosomal binding site in sequence.

        Arguments:
            sequence: the DNA sequence in upstream of the ORF to check for the RBS
            flexible: if true, searches for a flexible binding site (80% GA content within 8 bp), usefull for instnce in streptoycetes that oftentimes do not use a canonical rbs

        Returns:
            a boolean if the sequence contains a RBS
    """
    if RIBOSOMAL_BINDING_SITE.search(sequence) is not None:
        return True
    
    if flexible:
        if FLEXIBLE_RIBOSOMAL_BINDING_SITE.search(sequence) is not None:
            return True

    return False


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
        if cds.location.start + padding > last:
            intergenic_areas.append((max(start, last),
                                     min(end, int(cds.location.start) + padding)))
            last = int(cds.location.end) - padding
            continue
        # in case of existing CDS features overlapping, update the last end position
        if cds.location.start <= last <= cds.location.end:
            last = int(cds.location.end) - padding
    if last < end:
        intergenic_areas.append((max(start, last), end))
    return list(filter(lambda area: area[1] - area[0] >= min_length, intergenic_areas))


def find_all_orfs(record: Record, area: Optional[CDSCollection] = None,
                  min_length: int = 60, max_overlap: int = 10, include_rbs: bool = False, flexible_rbs: bool = False) -> List[CDSFeature]:
    """ Find ORFs within intergenic areas of the given record or subset of the record.

        Can (and should) be limited to just within a specific section of the record.

        Arguments:
            record: the record to search
            area: the specific CDSCollection to search within, or None
            min_length: the minimum length of ORFs to report, in nucleotides
            max_overlap: the maximum allowable bases of overlap with existing CDS features
            include_rbs: if true, search for an rbs and include that information in the feature

        Returns:
            a list of CDSFeatures, one for each ORF
    """
    # Get sequence for the range
    offset = 0
    seq = record.seq
    offset_end = len(record.seq)
    existing: Iterable[CDSFeature] = record.get_cds_features()
    if area:
        offset = area.location.start
        offset_end = area.location.end
        existing = record.get_cds_features_within_location(area.location,
                                                           with_overlapping=True)
    intergenic_areas = find_intergenic_areas(offset, offset_end, existing,
                                             min_length=min_length, padding=max_overlap)

    # Find orfs throughout the range
    locations = []
    for start, end in intergenic_areas:
        chunk = seq[start:end]
        locations.extend(scan_orfs(chunk, 1, start, minimum_length=min_length))
        locations.extend(scan_orfs(chunk.reverse_complement(), -1, start, minimum_length=min_length))
    new_features = []
    for location in locations:
        rbs = None
        if include_rbs:
            sequence = record.get_sequence_upstream_of_location(location, length=15)
            rbs = has_rbs(sequence, flexible=flexible_rbs)
        new_features.append(create_feature_from_location(record, location, rbs))

    return sorted(new_features)
