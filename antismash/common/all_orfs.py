# License: GNU Affero General Public License v3 or later
# A copy of GNU AGPL v3 should have been included in this software package in LICENSE.txt.

""" A very basic genefinding system. Finds all ORFs in a sequence, with some
    extra conditions applied.

    Not intended for genefinding an unannotated sequence, but for finding extra
    ORFs that have been skipped in an annotated sequence such as RiPP precursors.
"""


from typing import List, Optional
from typing import Iterable  # comment hints, pylint: disable=unused-import

from Bio.SeqFeature import FeatureLocation, BeforePosition, AfterPosition

from antismash.common.secmet import CDSFeature, Feature, Record
from antismash.common.secmet.features import CDSCollection


def scan_orfs(seq: str, direction: int, offset: int = 0) -> List[FeatureLocation]:
    """ Scan for open reading frames on a given sequence.
        Skips all ORFs with a size less than 60 bases.

        Arguments:
            seq: the sequence to examine
            direction: the search direction to use (all ORFs will use this as the strand)
            offset: an offset to add to any location discovered

        Returns:
            a list of FeatureLocations for each ORF, ordered by ascending position
    """
    seq = seq.upper()
    start_codons = ('ATG', 'GTG', 'TTG')
    stop_codons = ('TAA', 'TAG', 'TGA')
    matches = []
    # cache the sequence length
    seq_len = len(seq)
    for frame in [0, 1, 2]:
        i = frame
        last_stop = 0
        while i < seq_len - 2:
            if seq[i:i+3] in stop_codons and last_stop == 0:
                # special case for unstarted stops
                last_stop = i
                new_orf = FeatureLocation(BeforePosition(offset), offset + i + 2 + 1, direction)
                if direction == -1:
                    start = AfterPosition(seq_len + offset - new_orf.start)
                    end = seq_len + offset - new_orf.end
                    new_orf = FeatureLocation(end, start, strand=direction)
                matches.append(new_orf)
            if seq[i:i+3] not in start_codons:
                i += 3
                continue
            # Look for the next stop codon in this frame
            for j in range(i, seq_len - 2, 3):
                if seq[j:j+3] in stop_codons:
                    last_stop = j
                    # Skip Orfs that are shorter than 20 AA / 60 bases
                    if j - i <= 60:
                        break  # since no ORFs will be bigger before the stop
                    start = i
                    end = j + 2 + 1
                    if direction == 1:
                        new_orf = FeatureLocation(offset + start,
                                                  offset + end, direction)
                    else:
                        # reversed, so convert back to the forward positions
                        new_orf = FeatureLocation(seq_len + offset - end,
                                                  seq_len + offset - start, direction)
                    matches.append(new_orf)
                    # This was a good hit, update the last_stop cache.
                    break

            # if we found a matching stop, carry on looking for starts after this stop
            if last_stop > i:
                i = last_stop
                continue

            # Save orfs ending at the end of the sequence without stop codon
            if direction == 1:
                new_orf = FeatureLocation(i + offset, AfterPosition(seq_len + offset), direction)
            else:
                # reversed, so convert back to the forward positions
                new_orf = FeatureLocation(BeforePosition(offset), offset + seq_len - i, direction)
            matches.append(new_orf)
            # since there are no stop codons, just stop here
            break
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


def find_all_orfs(record: Record, area: Optional[CDSCollection] = None) -> List[CDSFeature]:
    """ Find all ORFs of at least 60 bases that don't overlap with existing
        CDS features.

        Can (and should) be limited to just within a specific section of the record.

        Arguments:
            record: the record to search
            area: the specific CDSCollection to search within, or None

        Returns:
            a list of CDSFeatures, one for each ORF
    """
    # Get sequence for the range
    offset = 0
    seq = record.seq
    existing = record.get_cds_features()  # type: Iterable[CDSFeature]
    if area:
        seq = area.extract(seq)
        offset = area.location.start
        existing = record.get_cds_features_within_location(area.location,
                                                           with_overlapping=True)

    # Find orfs throughout the range
    forward_matches = scan_orfs(seq, 1, offset)
    reverse_matches = scan_orfs(seq.reverse_complement(), -1, offset)
    locations = forward_matches + reverse_matches

    new_features = []

    for location in locations:
        if area:
            if isinstance(location.start, (BeforePosition, AfterPosition)):
                continue
            if isinstance(location.end, (BeforePosition, AfterPosition)):
                continue
        dummy_feature = Feature(location, feature_type="dummy")
        # skip if overlaps with existing CDSs
        if any(dummy_feature.overlaps_with(cds) for cds in existing):
            continue

        feature = create_feature_from_location(record, location)

        # skip if not wholly contained in the area
        if area and not feature.is_contained_by(area):
            continue

        new_features.append(feature)

    return new_features
