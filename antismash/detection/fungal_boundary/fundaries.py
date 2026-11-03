# License: GNU Affero General Public License v3 or later
# A copy of GNU AGPL v3 should have been included in this software package in
# LICENSE.txt.

""" Predict fungal BGC boundaries from promoter motif enrichment around core genes.

    Candidate boundaries must contain an antiSMASH core gene and a motif
    enriched relative to other motifs of its type on the same contig.
    Empirical significance is assessed against non-core genes outside the
    candidate windows (core gene neighbourhoob).

    Only predictions with an empirical p-value below config.P_VAL_MAX are
    reported. An empty report means that no prediction met the reporting
    requirements.
"""

from dataclasses import dataclass
from enum import IntEnum
from typing import Any

import numpy as np

from antismash.common.secmet import CDSFeature, GeneFunction, Record

from . import config_v7 as config
from . import utils_v7 as utils


@dataclass
class GeneWithMotifs:  # pylint: disable=too-many-instance-attributes
    """  One gene's annotation, promoter sequences and motif counts

        Coordinates are zero-based and end-exclusive. Both promoters follow
        the gene's reading direction. Each promoter-to-CDS distance measures
        bases from that promoter's start to the annotated 5' CDS boundary.
    """

    name: str
    start: int
    end: int
    strand: int
    promoter: str
    promoter_to_cds_start_dist: int
    linear_promoter: str
    linear_promoter_to_cds_start_dist: int
    motif_ids: np.ndarray
    is_anchor: bool
    counts: np.ndarray
    uncapped_counts: np.ndarray
    used_as_random_acting_anchor: bool = False

    @classmethod
    def from_feature(cls, record: Record, cds: CDSFeature) -> "GeneWithMotifs":
        """ Extracts CDS's promoters, identifies and counts motifs in the promoter

            Repeat motifs use up to config.UPSTREAM bases upstream of the
            CDS's 5' boundary.
            linear motifs use up to config.LINEAR_UPSTREAM
            of those bases. Both regions include up to config.DOWNSTREAM
            bases downstream and are limited by the contig boundaries.
            Promoters follow the gene's reading direction on either strand.
            Only CDS features assigned GeneFunction.CORE are anchors.

            Arguments:
                record: the record containing the CDS and DNA sequence of the contig
                cds: the CDS feature whose promoters are to be extracted

            Returns:
                the gene's annotation, promoters, and motif counts
        """
        cds_start_coord = int(cds.location.start)
        cds_end_coord = int(cds.location.end)
        strand = cds.location.strand

        if strand == -1:
            promoter_window_start = max(0, cds_end_coord - config.DOWNSTREAM)
            promoter_window_end = min(len(record.seq), cds_end_coord +
                                      config.UPSTREAM)
            promoter_seq = str(
                record.seq[promoter_window_start:promoter_window_end].
                reverse_complement()
            )
            promoter_to_cds_start_dist = promoter_window_end - cds_end_coord
        else:
            promoter_window_start = max(0, cds_start_coord - config.UPSTREAM)
            promoter_window_end = min(len(record.seq), cds_start_coord +
                                      config.DOWNSTREAM)
            promoter_seq = str(
                record.seq[promoter_window_start:promoter_window_end])
            promoter_to_cds_start_dist = (cds_start_coord
                                          - promoter_window_start)
        crop = max(0, promoter_to_cds_start_dist - config.LINEAR_UPSTREAM)
        linear_promoter = promoter_seq[crop:]
        motif_counts = count_motifs(promoter_seq, linear_promoter,
                                    capped_max_count=None)

        return cls(
            name=cds.get_name(),
            start=cds_start_coord,
            end=cds_end_coord,
            strand=strand,
            promoter=promoter_seq,
            promoter_to_cds_start_dist=promoter_to_cds_start_dist,
            linear_promoter=linear_promoter,
            linear_promoter_to_cds_start_dist=(
                promoter_to_cds_start_dist - crop),
            motif_ids=motif_counts.motif_ids,
            is_anchor=cds.gene_function == GeneFunction.CORE,
            counts=np.minimum(motif_counts.counts, config.MAX_COUNTS),
            uncapped_counts=motif_counts.counts,
        )


@dataclass
class ContigMotifs:
    """ Promoter evidence for one contig, indexed by gene name in CDS order.

        contig identifies the record. genes must preserve the record's CDS
        order, which defines gene adjacency and candidate boundary windows.
    """

    contig: str
    genes: dict[str, GeneWithMotifs]

    def get_gene(self, gene_name: str) -> GeneWithMotifs:
        """ Retrieves a gene's promoter information by its name.

            Arguments:
                gene_name: the name of a gene present on this contig

            Returns:
                the stored GeneWithMotifs instance
        """
        return self.genes[gene_name]


class SelectionMode(IntEnum):
    """ Motif proximity requirements for comparable enrichment scores.

        ANCHOR requires that the motif is in the anchor's promoter. 
        ADJACENT allows that a motif is found in immediate neighbour genes promoter.
        ALL allows a motif hit anywhere in the candidate window. 
    """

    ANCHOR = 0
    ADJACENT = 1
    ALL = 2


@dataclass(slots=True)
class EnrichedMotif:
    """ The most enriched motif found in a window of genes around one anchor.

        first_gene_offset and gene_count give the window relative to the
        anchor's own position, so it spans anchor_index + first_gene_offset
        through that plus gene_count. mode records how close to the anchor
        the motif was found: ANCHOR, ADJACENT, or ALL.
        p_empirical stays None until later where the calculate_empirical_p_values
        fills it in.
    """
    anchor_id: str
    motif_id: int
    llr: float
    mode: SelectionMode
    first_gene_offset: int
    gene_count: int
    p_empirical: float | None = None


@dataclass
class FundariesReport:
    """ Accepted boundaries and their supporting motif evidence for one contig.

        Each prediction contains anchor_id, start, end, motif_type,
        motif_string_dna, p_empirical, and gene_table. Each gene_table row
        contains name and motif_hits. All nested values are plain Python types.
        Coordinates are zero-based and end-exclusive.
    """

    contig: str
    predictions: list[dict[str, Any]]


def count_motifs(promoter: str, linear_promoter: str,
                 capped_max_count: int | None = 3) -> utils.MotifCounts:
    """Counts the repeat and linear motifs in a gene's two promoter.

        Repeat motifs are counted over the full promoter and linear motifs
        over the shorter linear promoter, then the two sets are concatenated.
        Repeat IDs come before linear IDs in the shared ID range, so the
        result stays sorted.

        Arguments:
            promoter: the full promoter sequence, used for repeat motifs
            linear_promoter: the shorter promoter used for linear motifs
            capped_max_count: the maximum count kept for each motif in
                              either promoter, or None to keep uncapped counts

        Returns:
            a MotifCounts instance covering both kinds of motif
    """
    repeat_counts = utils.repeat_motif_counts(
        promoter, capped_max_count=capped_max_count)
    linear_counts = utils.linear_motif_counts(
        linear_promoter, capped_max_count=capped_max_count)

    # Repeat IDs precede linear IDs, preserving sorted order here.
    return utils.MotifCounts(
        motif_ids=np.concatenate((repeat_counts.motif_ids,
                                  linear_counts.motif_ids)),
        counts=np.concatenate((repeat_counts.counts, linear_counts.counts)),
    )


def get_background_freqs_motif_ids(genes: dict[str, GeneWithMotifs]
                                   ) -> dict[int, float]:
    """ Returns each motif ID mapped to its frequency within its own type

        Counts are summed over every gene, using the capped counts, and
        each motif's total is then divided by the total for its motif type.
        The frequencies of each observed type therefore sum to 1, so a motif
        can only be enriched relative to other motifs of the same kind.

        Arguments:
            genes: the contig's genes in CDS order, keyed by gene name

        Returns:
            each observed motif ID mapped to its background frequency
    """
    background_freqs: dict[int, float] = {}
    occurrences_per_motif_type = {
        motif_type: 0.0 for motif_type in utils.MotifType
    }

    # Sum counts for each motif and each motif type across all genes.
    for gene in genes.values():
        for motif_id, count in zip(gene.motif_ids, gene.counts):
            motif_id = int(motif_id)
            count = float(count)
            motif_type = utils.get_motif_type(motif_id)
            background_freqs[motif_id] = (
                background_freqs.get(motif_id, 0.0) + count)
            occurrences_per_motif_type[motif_type] += count

    # Divide each motif's count by the total count for its type.
    for motif_id in background_freqs:
        motif_type = utils.get_motif_type(motif_id)
        background_freqs[motif_id] /= occurrences_per_motif_type[motif_type]

    return background_freqs


def find_the_most_enriched_motifs_within_window(
        genes, background_freqs_motif_ids, minimum_num_of_genes=3,
        max_genes_in_bgc=30, enriched_motif: EnrichedMotif | None = None,
) -> dict[str, EnrichedMotif | None]:
    """ Finds the most enriched motif in a window of genes around each anchor

        Runs in one of two modes. With enriched_motif left as None, every
        motif is scored at every gene antiSMASH called core gene. Given an
        enriched_motif, only that motif is scored, only at the genes marked
        as random acting anchors, and only in its own mode, which builds the
        null distribution for an earlier result.

        Windows are ranked by a Poisson log-likelihood ratio of observed
        against expected counts, the expectation being the motif's background
        frequency times the total count of its motif type in the window. To be
        eligible a motif must occur in at least minimum_num_of_genes genes of
        the window and exceed its expected count. Repeat motifs are preferred
        over linear ones only when their LLR scores are equal. Remaining ties
        are resolved by mode, then the fewest genes in the window, then the
        lowest motif ID.

        Arguments:
            genes: the contig's genes in CDS order, keyed by gene name
            background_freqs_motif_ids: each motif ID mapped to its background
                                        frequency within its own motif type
            minimum_num_of_genes: the fewest genes of a window a motif must
                                 occur in to be considered
            max_genes_in_bgc: the most genes a window may span, and the
                              furthest it may reach from its anchor
            enriched_motif: an earlier result to rescore for a null
                            distribution, or None to search all motifs
                            at the real anchors

        Returns:
            each acting anchor's name mapped to its best EnrichedMotif, or to
            None where no motif was eligible
    """
    if enriched_motif is None:
        # Search all possible motif IDs at the real anchors.
        # len(motif_ids) >= 1
        motif_ids = np.arange(config.N_MOTIFS)
    else:
        # Find the most enriched window for the supplied motif.
        # len(motif_ids) == 1
        motif_ids = np.array([enriched_motif.motif_id], dtype=int)

    # To find out what kind of motif_type (inverted, direct, or linear)
    # the motif_ids is 
    motif_types = np.searchsorted(config.MOTIF_TYPE_ID_BOUNDARIES, motif_ids,
                                  side="right")

    if enriched_motif is None:
        # Search all promoters before classifying the winning motif's mode.
        required_mode = SelectionMode.ALL
    else:
        required_mode = enriched_motif.mode
    
    gene_list = list(genes.values())
    results = {}
    for anchor_index, anchor in enumerate(gene_list):
        # treat "used_as_random_acting_anchor" and "is_anchor" as acting_anchor
        acting_anchor = (anchor.is_anchor if enriched_motif is None
                         else anchor.used_as_random_acting_anchor)
        if not acting_anchor:
            continue
        
        results[anchor.name] = None
        if enriched_motif is not None and required_mode < SelectionMode.ALL:
            if required_mode == SelectionMode.ANCHOR:
                # The motif must occur in the anchor itself.
                genes_to_check = [anchor]
            else:
                # The motif may occur in the anchor or either adjacent gene.
                first_index = max(0, anchor_index - 1)
                stop_index = anchor_index + 2
                genes_to_check = gene_list[first_index:stop_index]
            # at least one gene in genes_to_check contains the target motif
            motif_found = any(
                enriched_motif.motif_id in gene.motif_ids
                for gene in genes_to_check
            )
            if not motif_found:
                continue
        # making sure that the window is inside the contig/gene_list
        start_gene_index_in_window = max(0, anchor_index - max_genes_in_bgc)
        end_gene_index_in_window = min(len(gene_list), anchor_index + 
                                       max_genes_in_bgc + 1)
        anchor_index_in_window = anchor_index - start_gene_index_in_window
        # initalize windows lenght x num_motifs 
        # which stores the count data, for the genes in the current window
        # that is gene_0 has 2 occurances of motif_id 10
        counts = np.zeros((end_gene_index_in_window -
                           start_gene_index_in_window,
                           len(motif_ids)), dtype=np.int64)
        # initalize windows_lenght x 3 (different kind of motif types)
        type_counts = np.zeros((len(counts), 3), dtype=np.int64)
        # Loop through each gene object and its zero-based position within
        # the window.
        for gene_index_in_window, gene in enumerate(
                gene_list[start_gene_index_in_window:end_gene_index_in_window]):
            if enriched_motif is None:
                counts[gene_index_in_window, gene.motif_ids] = gene.counts
            else:
                counts[gene_index_in_window, 0] = (gene.counts[gene.motif_ids == enriched_motif.motif_id].sum())
            # To find out what kind of motif_type (inverted, direct, or linear) the motif_ids is 
            gene_types = np.searchsorted(config.MOTIF_TYPE_ID_BOUNDARIES, 
                                         gene.motif_ids, side="right")
            # adds up motif-id occurrences by motif-type for the current gene
            # np.bincount: Count number of occurrences of each value in array
            type_counts[gene_index_in_window] = np.bincount(
                gene_types, weights=gene.counts, minlength=3)

        if (enriched_motif is not None
                and np.count_nonzero(counts) < minimum_num_of_genes):
            continue
        # initialize a 2d matrix: (number of genes in the window +1) x motif IDs
        # +1 extra row is used for cumulative sums later
        cumulative_occurrences_per_motif = np.zeros(
            (len(counts) + 1, len(motif_ids)), dtype=np.int64)
        cumulative_genes_with_motif = (
            np.zeros_like(cumulative_occurrences_per_motif))
        # Store cumulative motif occurrence counts by type:
        # inverted, direct, linear.
        cumulative_occurrences_per_motif_type = np.zeros(
            (len(counts) + 1, 3), dtype=np.int64)
        # Calculate running totals for each motif across genes.
        # Example, using these cols = motif IDs 10, 42, 80
        # counts       cumulative_occurrences_per_motif
        #              [0, 0, 0]  <- no genes counted yet
        # [2, 0, 1]    [2, 0, 1]  <- gene 0
        # [1, 3, 0]    [3, 3, 1]  <- genes 0 and 1
        # [0, 1, 2]    [3, 4, 3]  <- all three genes
        np.cumsum(counts, axis=0, out=cumulative_occurrences_per_motif[1:])
        np.cumsum(counts > 0, axis=0, out=cumulative_genes_with_motif[1:])
        np.cumsum(type_counts, axis=0, 
                  out=cumulative_occurrences_per_motif_type[1:])
        in_anchor = counts[anchor_index_in_window] > 0

        best_rank = None
        # Loop which tries different starting genes for candidate
        # windows containing the anchor gen
        # in most cases (not right next to contig end), will the range be from 0 to 31 
        for first in range(
            max(0, anchor_index_in_window - max_genes_in_bgc + 1),
            anchor_index_in_window + 1
        ):
            # first_stop and last_stop defines the earliest and latest allowed index for window
            first_stop = max(anchor_index_in_window + 1,
                             first + minimum_num_of_genes)
            last_stop = min(len(counts), first + max_genes_in_bgc)
            stops = np.arange(first_stop, last_stop + 1)

            # the idea is to look at the 
            # Total inside a candidate window = running total at its end − running total before its start.
            observed = (cumulative_occurrences_per_motif[stops]
                        - cumulative_occurrences_per_motif[first])
            genes_with_motif = (cumulative_genes_with_motif[stops]
                                - cumulative_genes_with_motif[first])
            window_type_totals = (cumulative_occurrences_per_motif_type[stops]
                                  - cumulative_occurrences_per_motif_type[first])
            # expected count = total occurrences of its motif type in the window
            #     × its background frequency within that type
            # Expected occurrences of each motif in each candidate window:
            # multiply its motif-type total by its background frequency within that type.
            expected = window_type_totals[:, motif_types] * np.array([
                background_freqs_motif_ids.get(motif_id, 0.0)
                for motif_id in motif_ids
            ], dtype=float)
            nearby_stops = np.minimum(stops, anchor_index_in_window + 2)
            in_adjacent = (
                cumulative_genes_with_motif[nearby_stops]
                - cumulative_genes_with_motif[max(first,
                                                  anchor_index_in_window - 1)]
            ) > 0

            eligible = ((genes_with_motif >= minimum_num_of_genes) &
                        (observed > expected))
            if required_mode == SelectionMode.ANCHOR:
                eligible &= in_anchor
            elif required_mode == SelectionMode.ADJACENT:
                eligible &= in_adjacent
            rows, columns = np.nonzero(eligible)
            if columns.size == 0:
                continue
            hits = observed[rows, columns]
            expected_hits = expected[rows, columns]
            scores = (hits * np.log(hits / expected_hits) -
                      (hits - expected_hits))
            highest_llr = float(scores.max())
            for winner in np.flatnonzero(scores == highest_llr):
                row, column = rows[winner], columns[winner]
                stop = int(stops[row])
                motif_id = int(motif_ids[column])
                start = first - anchor_index_in_window
                gene_count = stop - first
                if enriched_motif is not None:
                    mode_rank = required_mode
                elif in_anchor[column]:
                    mode_rank = SelectionMode.ANCHOR
                elif in_adjacent[row, column]:
                    mode_rank = SelectionMode.ADJACENT
                else:
                    mode_rank = SelectionMode.ALL
                is_linear = motif_id >= config.LINEAR_MOTIF_ID_START

                # Compare priorities if different candidatas achieve
                # the same score.
                candidate_rank = (
                    -highest_llr,   # Higher LLR first.
                    is_linear,      # Repeats (False) win when LLR scores tie.
                    mode_rank,      # Anchor, then adjacent, then all.
                    gene_count,     # Fewer genes in the window first.
                    motif_id,       # Lower motif ID breaks the final tie.
                )
                # Keep the first candidate when all ranking values are equal.
                if best_rank is None or candidate_rank < best_rank:
                    best_rank = candidate_rank
                    results[anchor.name] = EnrichedMotif(
                        anchor_id=anchor.name,
                        motif_id=motif_id,
                        llr=highest_llr,
                        mode=mode_rank,
                        first_gene_offset=start,
                        gene_count=gene_count,
                    )
    return results


def _eligible_anchor_positions(gene_list: list[GeneWithMotifs],
                               enriched_motifs: dict[str, EnrichedMotif | None],
                               reach: int) -> np.ndarray:
    """ Identifies non-core controls whose spans avoid candidate cluster windows.

        A gene is eligible when antiSMASH did not call it a core gene, it has
        max_genes_in_bgc genes either side of it within the contig, and no
        gene in that span is already covered by one of the windows in
        enriched_motifs.

        Arguments:
            gene_list: the contig's genes in CDS order
            enriched_motifs: candidate results whose windows exclude control
                             genes, regardless of score or p-value
            reach: the required number of neighbouring genes on each side
                   that must lie within the contig and outside candidate windows

        Returns:
            an array of positions into gene_list, in increasing order
    """
    gene_indices = {gene.name: index for index, gene in enumerate(gene_list)}

    # A gene cannot act as an anchor if antiSMASH already called it one, or if
    # its own window would reach a gene that a reported window sits on. The
    # second test is the same as the gene lying within reach of such a window,
    # so widening each reported window by reach blocks exactly those genes.
    is_blocked = np.array([gene.is_anchor for gene in gene_list], dtype=bool)
    for result in enriched_motifs.values():
        if result is None:
            continue
        # first_gene_offset counts from the anchor itself and is negative when
        # the window opens before it, so adding it gives an absolute position.
        window_start = (gene_indices[result.anchor_id] +
                        result.first_gene_offset)
        window_end = window_start + result.gene_count
        # max() because a negative slice bound would wrap to the contig's end
        is_blocked[max(0, window_start - reach):window_end + reach] = True

    # Only a gene with a full window's worth of neighbours on either side can
    # act as an anchor, which rules out the first and last few on the contig.
    positions = np.arange(reach, len(gene_list) - reach)
    return positions[~is_blocked[positions]]


def _sample_spaced_positions(candidates: np.ndarray, spacing:
                             int, max_count: int, rng: np.random.Generator
                             ) -> np.ndarray:
    """ Draws a random subset of candidates that are spaced apart

        Candidates are visited in a random order and kept when they clear
        every candidate kept so far, so the spacing always holds. Fewer than
        max_count are returned when the spacing does not allow that many,
        which on a real contig is the usual case.

        The result is a maximal selection rather than the largest possible
        one: an early pick can block two later ones that would have fitted
        together, so a luckier order would sometimes fit more.

        Arguments:
            candidates: the sorted positions to choose from
            spacing: the fewest positions two chosen candidates may differ by
            max_count: the most candidates to choose
            rng: the generator to draw from

        Returns:
            an array of the chosen positions, in increasing order
    """
    # Basically chose different anchors that are separated in `spacing´
    # so we are selecting valuable anchors and not too close ones.
    chosen: list[int] = []
    for index in rng.permutation(len(candidates)):
        if len(chosen) >= max_count:
            break
        position = int(candidates[index])
        if all(abs(position - taken) >= spacing for taken in chosen):
            chosen.append(position)
    return np.array(sorted(chosen), dtype=np.int64)


def select_random_acting_anchors(
        genes: dict[str, GeneWithMotifs],
        enriched_motifs: dict[str, EnrichedMotif | None],
        max_genes_in_bgc: int = 30) -> None:
    """ Marks a spaced, random sample of genes to act as anchors for the
        calculation of p_empirical.

        Sets used_as_random_acting_anchor in place, clearing any previous
        selection first. Candidates come from _eligible_anchor_positions and
        are spaced so that no two of their windows overlap by more than
        config.MAX_RANDOM_ANCHOR_OVERLAP genes. At most
        config.MAX_RANDOM_ANCHORS are marked, and fewer when the spacing does
        not allow that many, which is the usual case: how many fit is set by
        the spacing rather than by that ceiling.

        The sample is drawn by visiting the candidates in a random order
        fixed by config.RANDOM_SEED and keeping each one that is far enough
        from those already kept.

        Arguments:
            genes: the contig's genes in CDS order
            enriched_motifs: each anchor's result from
                             find_the_most_enriched_motifs_within_window,
                             whose windows are excluded from the candidates
            max_genes_in_bgc: the most genes a window may span

        Returns:
            None; the selection is recorded on the genes themselves
    """
    gene_list = list(genes.values())
    for gene in gene_list:
        gene.used_as_random_acting_anchor = False

    # A window reaches this far either side of its anchor, so two windows
    # sharing at most config.MAX_RANDOM_ANCHOR_OVERLAP genes must have their
    # anchors at least this far apart.
    reach = max_genes_in_bgc - 1
    # config.MAX_RANDOM_ANCHOR_OVERLAP and reach should maybe be the same?
    spacing = max(1, 2 * reach + 1 - config.MAX_RANDOM_ANCHOR_OVERLAP)
    candidates = _eligible_anchor_positions(gene_list, enriched_motifs, reach)

    chosen = _sample_spaced_positions(candidates, spacing,
                                      config.MAX_RANDOM_ANCHORS,
                                      np.random.default_rng(
                                          config.RANDOM_SEED))
    for position in chosen:
        gene_list[int(position)].used_as_random_acting_anchor = True


def calculate_empirical_p_values(
        genes: dict[str, GeneWithMotifs],
        enriched_motifs: dict[str, EnrichedMotif | None],
        background_freqs_motif_ids: dict[int, float]) -> None:
    """ Fills in an empirical p-value for each enriched motif result

        p_empirical is filled in only for the results that are tested, and
        left as None on the rest.Only results with an LLR of at least
        config.MIN_LLR are tested. For each distinct motif and mode among
        those, that motif is rescored at the random acting anchors to build
        a null distribution, and the p-value is the share of those null
        scores reaching the observed LLR.
        One is added to both numerator and denominator, so a p-value is
        never exactly zero..

        Nothing is computed when fewer than config.MIN_RANDOM_ANCHORS random
        acting anchors were selected, and a warning is flagged whenever there
        are fewer than config.MAX_RANDOM_ANCHORS of them.

        Arguments:
            genes: the contig's genes in CDS order, keyed by gene name
            enriched_motifs: the results to annotate, as returned by
                             find_the_most_enriched_motifs_within_window
            background_freqs_motif_ids: each motif ID mapped to its background
                                        frequency within its own motif type

        Returns:
            None; p_empirical is set on the results themselves
    """
    for result in enriched_motifs.values():
        # Guard for the anchors that did not have any alid motifs
        if result is not None:
            result.p_empirical = None
    sample_size = sum(gene.used_as_random_acting_anchor
                      for gene in genes.values())
    if sample_size < config.MAX_RANDOM_ANCHORS:
        print(f"Warning: only {sample_size} random acting anchors available "
              f"(maximum {config.MAX_RANDOM_ANCHORS}).")
    if sample_size < config.MIN_RANDOM_ANCHORS:
        print("Skipping empirical p-values: fewer than "
              f"{config.MIN_RANDOM_ANCHORS} random acting anchors.")
        return

    random_scores = {}
    for result in enriched_motifs.values():
        if result is None or result.llr < config.MIN_LLR:
            continue
        key = (result.motif_id, result.mode)
        if key not in random_scores:
            random_results = find_the_most_enriched_motifs_within_window(
                genes, background_freqs_motif_ids, enriched_motif=result
            )
            random_scores[key] = np.array([
                0.0 if random_result is None else random_result.llr
                for random_result in random_results.values()
            ])
        result.p_empirical = float(
            (1 + np.count_nonzero(random_scores[key] >= result.llr))
            / (sample_size + 1))


def build_report(contig: str, genes: dict[str, GeneWithMotifs],
                 enriched_motifs: dict[str, EnrichedMotif | None]
                 ) -> FundariesReport:
    """ Collects the accepted results into a report for one contig.

        Only results with an empirical p-value below config.P_VAL_MAX are
        reported. Each accepted window is located from its anchor's
        position plus the result's first_gene_offset, spanning gene_count
        genes. Every gene row stores the uncapped count.

        Arguments:
            contig: the record ID to store on the report
            genes: the contig's genes in CDS order, keyed by gene name
            enriched_motifs: each anchor's result

        Returns:
            a FundariesReport for this contig, holding one prediction per
    """
    gene_list = list(genes.values())
    gene_indices = {gene.name: index for index, gene in enumerate(gene_list)}
    predictions: list[dict[str, Any]] = []
    for result in enriched_motifs.values():
        if (result is None or result.p_empirical is None
                or result.p_empirical >= config.P_VAL_MAX):
            continue

        first_index = gene_indices[result.anchor_id] + result.first_gene_offset
        window_genes = gene_list[first_index:first_index + result.gene_count]
        motif_type = utils.get_motif_type(result.motif_id)
        gene_table = []
        for gene in window_genes:
            motif_hits = int(gene.uncapped_counts[
                gene.motif_ids == result.motif_id].sum())
            gene_table.append({"name": gene.name, "motif_hits": motif_hits})

        predictions.append({
            "anchor_id": result.anchor_id,
            "start": min(gene.start for gene in window_genes),
            "end": max(gene.end for gene in window_genes),
            "motif_type": motif_type.value,
            "motif_string_dna": utils.motif_int_to_text(result.motif_id),
            "p_empirical": result.p_empirical,
            "gene_table": gene_table,
        })
    return FundariesReport(contig=contig, predictions=predictions)


def main(record: Record) -> FundariesReport:
    """ Reads a record from antiSMASH, scores every anchor and returns
        the accepted boundaries.

        Each CDS becomes a GeneWithMotifs, then every gene antiSMASH
        called an core is scored for the most enrichsd motif. A sample
        of other genes is used to re-score the same motif to build an
        empirical p-value. Only the motifs beating config.P_VAL_MAX 
        reach the report.

        Arguments:
            record: the record whose CDS features will be scored

        Returns:
            a FundariesReport containing the contig ID and predictions with
            available empirical p-values below config.P_VAL_MAX, or an empty
            report if the contig has fewer than config.MIN_GENES_PER_CONTIG
            CDS features
    """
    cds_features = record.get_cds_features()
    if len(cds_features) < config.MIN_GENES_PER_CONTIG:
        return FundariesReport(contig=record.id, predictions=[])

    genes = {}
    for cds in cds_features:
        gene = GeneWithMotifs.from_feature(record, cds)
        genes[gene.name] = gene

    background_freqs_motif_ids = get_background_freqs_motif_ids(genes)
    results = find_the_most_enriched_motifs_within_window(
        genes=genes,
        background_freqs_motif_ids=background_freqs_motif_ids,
                )
    
    select_random_acting_anchors(genes, results)
    calculate_empirical_p_values(genes, results, background_freqs_motif_ids)

    return build_report(record.id, genes, results)
