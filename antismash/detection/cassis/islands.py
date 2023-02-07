# License: GNU Affero General Public License v3 or later
# A copy of GNU AGPL v3 should have been included in this software package in LICENSE.txt.

""" Island specific logic and classes for CASSIS """

import logging
from typing import Any, List

from .config import MAX_GAP_LENGTH, VERBOSE_DEBUG
from .motifs import Motif
from .promoters import Promoter


class Island:
    """ A container for an island between two promoters with a motif. """
    def __init__(self, start: Promoter, end: Promoter, motif: Motif) -> None:
        assert isinstance(start, Promoter)
        self.start = start
        assert isinstance(end, Promoter)
        self.end = end
        assert isinstance(motif, Motif)
        self.motif = motif

    def __eq__(self, other: Any) -> bool:
        return (isinstance(other, Island)
                and self.start == other.start
                and self.end == other.end
                and self.motif == other.motif)

    def __repr__(self) -> str:
        return f"Island(start={self.start}, end={self.end}, motif={self.motif})"


def get_islands(anchor_promoter: int, motifs: List[Motif], promoters: List[Promoter]) -> List[Island]:
    """Find islands of binding sites (previously found by FIMO) around anchor gene to define cluster borders"""
    assert MAX_GAP_LENGTH >= 0
    islands = []
    motifs = list(motifs)
    for motif in motifs:
        # create list with binding sites per promoter
        bs_per_promoter = [0] * len(promoters)  # first: set number of binding sites to 0
        for i, promoter in enumerate(promoters):
            # second: set actual number of binding sites, if any
            bs_per_promoter[i] = motif.hits[promoter.get_id()]

        # upstream
        start = anchor_promoter  # init upstream cluster border
        i = anchor_promoter  # init position of anchor gene's promoter
        while i > 0:
            # promoter with binding site
            # … 1 …
            if bs_per_promoter[i-1] > 0:
                start -= 1
                i -= 1
                continue
            # check with gaps, e.g. 1 0 1 or 1 0 0 1
            gap_size = get_island_gap_size(i, bs_per_promoter, MAX_GAP_LENGTH, upstream=True)
            if gap_size == 0:
                break
            start -= gap_size + 1
            i -= gap_size + 1

        # downstream
        i = anchor_promoter  # reset position of anchor gene's promoter
        end = anchor_promoter  # init downstream cluster border
        while i < len(bs_per_promoter) - 1:

            # promoter with binding site(s)
            # … 1 …
            if bs_per_promoter[i + 1] > 0:
                end += 1
                i += 1
                continue

            # check with gaps, e.g. 1 0 1 or 1 0 0 1
            gap_size = get_island_gap_size(i, bs_per_promoter, MAX_GAP_LENGTH, upstream=False)
            if gap_size == 0:
                break
            end += gap_size + 1
            i += gap_size + 1

        if VERBOSE_DEBUG:
            logging.debug("Island %s -- %s (motif %s)", promoters[start].get_id(),
                          promoters[end].get_id(), motif)
        islands.append(Island(promoters[start], promoters[end], motif))
    return islands


def get_island_gap_size(position: int, promoter_scores: List[int], max_gap_size: int,
                        upstream: bool) -> int:
    """ Finds the size of an island in promoters, if it exists.

        Arguments:
            position: the position around which to check for an island in promoter_scores
            promoter_scores: the number of hits for each promoter, ordered by location
            max_gap_size: the maximum gap size allowable for an island
            upstream: whether the search is looking upstream or not

        Returns:
            0 if no island found, otherwise an int less than or equal to max_gap_size
    """
    if upstream:
        end = position + 1
    else:
        start = position

    for gap_size in range(1, max_gap_size + 1):
        if upstream:
            start = position - gap_size - 1
        else:
            end = position + gap_size + 2

        # skip anything that there are no promoters for
        if start < 0 or end > len(promoter_scores):
            continue

        # for gap_size of 1: x 0 x, where 0 is the anchor
        # for gap_size of 2: x 0 0 x...
        # and so on
        chunk = promoter_scores[start:end]
        if chunk[0] > 0 and chunk[-1] > 0 and chunk[1:-1] == [0] * gap_size:
            # since no larger island can contain a smaller one, stop
            return gap_size
    return 0
