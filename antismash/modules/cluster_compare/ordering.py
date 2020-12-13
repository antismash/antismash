# License: GNU Affero General Public License v3 or later
# A copy of GNU AGPL v3 should have been included in this software package in LICENSE.txt.

""" Contains functions for calculating the ordering score of a reference/query pairing.

The order score takes into account how much of the pairing is contiguous in
both query and reference. Scoring penalties apply for gaps and for strand
mismatches of chunks.
"""


from typing import (
    Dict,
    List,
    Sequence,
    Tuple,
)

from antismash.common.secmet import CDSFeature

from .data_structures import Hit, ReferenceCDS

Pairing = Tuple[int, int, bool]


NORMALISATION_TARGET = (0.8, 0.8)  # (desired score, required percentage contiguous matches)
EXTRA_SEGMENT_PENALTY = 0.75  # for multiple disjoint segments
REVERSED_SEGMENT_PENALTY = 0.9  # for when the strand of a segment doesn't match


def calculate_order_score(area_features: Sequence[CDSFeature], hits: Dict[str, Hit],
                          references: Dict[str, ReferenceCDS]) -> float:
    """ Calculates a score (0..1) of the relative ordering of linked features within
        the query area and the reference area. A score of 1 requires that all
        features in one of the two areas have links and that the order of the
        linked features is consistent.

        Arguments:
            area_features: a sequence of CDSFeatures from the query area
            hits: a dictionary mapping reference CDS name to matching Hit instance
            references: a dictionary mapping reference CDS name to ReferenceCDS instance

        Returns:
            a float ranging from 0 to 1, inclusive
    """
    if not hits:
        return 0.

    segments = find_segments(hits, area_features, references)
    assert segments

    return score_segments(segments, len(area_features), len(references))


def _build_segments_from_pairings(pairings: Sequence[Pairing]) -> List[List[Pairing]]:
    """ Takes a sequence of pairings and separates them into contiguous chunks.

        A mismatching strand is considered a break in contiguity of matching strands.

        Arguments:
            pairings: a sequence of tuples of the form (query index, reference index, strand match)

        Returns:
            all provided pairings in a series of lists, each list being a contiguous chunk
    """
    if not pairings:
        return []
    segments = [[pairings[0]]]
    for pairing in pairings[1:]:
        prev_cds, prev_ref, prev_strand_match = segments[-1][-1]
        cds, ref, strand_match = pairing
        if cds <= prev_cds:
            raise ValueError("pairs not ordered by first index")
        if any([
            prev_strand_match != strand_match,  # strands no longer match
            cds != prev_cds + 1,  # CDSes aren't contiguous
            abs(ref - prev_ref) != 1,   # references aren't contiguous
        ]):
            segments.append([])
        segments[-1].append(pairing)
    return segments


def find_segments(hits: Dict[str, Hit], features: Sequence[CDSFeature],
                  reference_features: Dict[str, ReferenceCDS]) -> List[List[Pairing]]:
    """ Finds segments of contiguous hits from the given hits and features

        Arguments:
            hits: a dictionary mapping reference ID to Hit instance
            features: a list of CDSFeature instance, sorted by location
            reference_features: a dictionary mapping reference ID to ReferenceCDS instance

        Returns:
            a series of lists, each list being a contiguous chunk of hits, with
            each hit a tuple in the form (query index, reference index, strand match)
    """
    pairings = []
    for ref_cds_name, hit in hits.items():
        cds_index = features.index(hit.cds) + 1
        reference_index = int(hit.reference_id)
        reference_strand = reference_features[ref_cds_name].location.strand
        pairings.append((cds_index, reference_index, hit.cds.location.strand == reference_strand))

    pairings.sort()

    return _build_segments_from_pairings(pairings)


def score_segments(segments: List[List[Pairing]], query_size: int, reference_size: int) -> float:
    """ Calculates the normalised score of the given segments. Normalisation
        is based on the theoretical best case.

        Arguments:
            segments: a series of lists, each list being a contiguous chunk of hits
            query_size: the number of CDS features in the query
            reference_size: the number of CDS features in the reference

        Returns:
            a float ranging between 0 and 1, inclusive
    """
    max_possible = min(query_size, reference_size)

    # calculate the base of the exponential scoring function
    base = (1/NORMALISATION_TARGET[1])**(1/((1-NORMALISATION_TARGET[0])*max_possible))

    result = base**sum(len(segment) for segment in segments)

    # penalise for each extra segment based on the gap size (in genes)
    for i, segment in enumerate(segments[1:]):
        # find the distance between each segment for the query
        query_gap = (segment[0][0] - segments[i][-1][0]) / query_size
        # account for references not being ordered
        if segments[i][0][1] < segment[0][1]:  # increasing order
            prev_ref = max(segments[i][0][1], segments[i][-1][1])
            current_ref = min(segment[0][1], segment[-1][1])
            ref_gap = (current_ref - prev_ref) / reference_size
        else:  # decreasing order
            prev_ref = min(segments[i][0][1], segments[i][-1][1])
            current_ref = max(segment[0][1], segment[-1][1])
            ref_gap = (prev_ref - current_ref) / reference_size
        # then use the absolute difference in case they're reversed
        result *= 0.5 + 0.5 * (1 - max(query_gap, ref_gap))

    # penalise for having incorrect strands for a segment
    reversed_segments = 0
    for segment in segments:
        if not segment[0][2]:  # strand mismatch
            reversed_segments += 1
    # since the whole query could be reversed, take the complement if relevant
    reversed_segments = min(reversed_segments, len(segments) - reversed_segments)
    result *= REVERSED_SEGMENT_PENALTY**reversed_segments

    # normalise by theoretical perfect score
    result = result / base**max_possible
    return result
