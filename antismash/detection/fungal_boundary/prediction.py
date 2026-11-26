# License: GNU Affero General Public License v3 or later
# A copy of GNU AGPL v3 should have been included in this software package in
# LICENSE.txt.

""" Predict fungal BGC boundaries from promoter motif enrichment around core genes.

    Candidate boundaries must contain a core gene and a motif
    enriched relative to other motifs of its type on the same contig.
    Boundaries are reported when motif enrichment meets the minimum
    log-likelihood ratio (LLR).
"""

from __future__ import annotations

from collections.abc import Iterator
from dataclasses import asdict, dataclass
from enum import IntEnum
from typing import Any, Self

import numpy as np

from antismash.common.module_results import DetectionResults
from antismash.common.secmet import CDSFeature, FeatureLocation, GeneFunction, Record, SubRegion

from . import utils
from .config import config


@dataclass
class GeneWithMotifs:  # pylint: disable=too-many-instance-attributes
    """ One gene's basic annotations like promoter sequences and motif counts
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
    capped_motif_counts: np.ndarray  # capped per motif within a single promoter
    uncapped_motif_counts: np.ndarray

    @classmethod
    def from_feature(cls, record: Record, cds: CDSFeature) -> GeneWithMotifs:
        """ Extracts the CDS’s promoters, identifies motifs, and counts occurrences
            of each motif.

            Repeat and linear motif promoters extend up to config.UPSTREAM_BASES
            and config.LINEAR_UPSTREAM_BASES bases upstream of the CDS's 5' boundary,
            respectively. Both extend up to config.DOWNSTREAM_BASES bases downstream,
            stop at contig boundaries, and follow the gene's reading direction.
            Only CDS features assigned GeneFunction.CORE are anchors.

            Arguments:
                record: the record containing the CDS and DNA sequence of the contig
                cds: the CDS feature whose promoters are to be extracted

            Returns:
                the gene's annotations
        """
        cds_start_coord = int(cds.location.start)
        cds_end_coord = int(cds.location.end)
        strand = cds.location.strand

        if strand == -1:
            promoter_window_start = max(0, cds_end_coord - config.DOWNSTREAM_BASES)
            promoter_window_end = min(len(record.seq), cds_end_coord +
                                      config.UPSTREAM_BASES)
            promoter_seq = str(
                record.seq[promoter_window_start:promoter_window_end].
                reverse_complement()
            )
            promoter_to_cds_start_dist = promoter_window_end - cds_end_coord
        else:
            promoter_window_start = max(0, cds_start_coord - config.UPSTREAM_BASES)
            promoter_window_end = min(len(record.seq), cds_start_coord +
                                      config.DOWNSTREAM_BASES)
            promoter_seq = str(
                record.seq[promoter_window_start:promoter_window_end])
            promoter_to_cds_start_dist = (cds_start_coord
                                          - promoter_window_start)
        crop = max(0, promoter_to_cds_start_dist - config.LINEAR_UPSTREAM_BASES)
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
            capped_motif_counts=np.minimum(
                motif_counts.counts, config.MAX_MOTIFS_PER_PROMOTER),
            uncapped_motif_counts=motif_counts.counts,
        )


class SelectionMode(IntEnum):
    """ Defines the modes - or where a motif must occur relative to the anchor gene.

        ANCHOR requires that the motif is in the anchor's promoter.
        ADJACENT allows the motif to occur in the promoter of either
        immediate neighbour of the anchor gene.
        ALL allows a motif hit anywhere in the candidate window.
    """

    ANCHOR = 0
    ADJACENT = 1
    ALL = 2


@dataclass(slots=True)
class EnrichedMotif:
    """ Finding the most enriched motif and its gene window around an anchor.

        first_gene_offset gives the window's start relative to the anchor:
        negative values count genes before it, positive values count genes after it,
        and zero means the window starts at the anchor.
        gene_count is the number of genes in the window.
        mode describes the motif's proximity to the anchor gene.
    """
    anchor_id: str
    motif_id: int
    llr: float
    mode: SelectionMode
    first_gene_offset: int
    gene_count: int


@dataclass
class BoundaryPrediction:
    """ An accepted boundary and its supporting motif evidence.

        motif_counts maps gene names to uncapped/true motif counts in CDS order.
    """
    anchor_id: str
    start: int
    end: int
    motif_type: str
    motif_string_dna: str
    motif_counts: dict[str, int]

    def to_json(self) -> dict[str, Any]:
        """Serialise a BoundaryPrediction"""
        return asdict(self)

    @classmethod
    def from_json(cls, json: dict[str, Any]) -> Self:
        return cls(
            anchor_id=json["anchor_id"],
            start=json["start"],
            end=json["end"],
            motif_type=json["motif_type"],
            motif_string_dna=json["motif_string_dna"],
            motif_counts=dict(json["motif_counts"]),
        )


@dataclass
class BoundaryDetectionResults(DetectionResults):
    """ Accepted boundaries and their supporting motif evidence for one contig """

    schema_version = 2
    contig: str
    predictions: list[BoundaryPrediction]

    def __init__(self, contig: str, predictions: list[BoundaryPrediction]) -> None:
        super().__init__(contig)
        self.contig = contig
        self.predictions = predictions

    def to_json(self) -> dict[str, Any]:
        """Serialize the report and its schema version"""
        return {
            "schema_version": self.schema_version,
            "record_id": self.record_id,
            "report": {
                "contig": self.contig,
                "predictions": [p.to_json() for p in self.predictions],
            },
        }

    @staticmethod
    def from_json(json: dict[str, Any], record: Record) -> BoundaryDetectionResults | None:
        """Reconstruct boundary detection results for the original record."""
        if json.get("schema_version") != BoundaryDetectionResults.schema_version:
            return None
        if json.get("record_id") != record.id:
            return None
        report_data = json["report"]
        predictions = [
            BoundaryPrediction.from_json(data) for data in report_data.get("predictions", [])
        ]
        return BoundaryDetectionResults(
            contig=report_data["contig"],
            predictions=predictions,
        )

    def get_predicted_subregions(self) -> list[SubRegion]:
        """Create subregions for the predicted boundaries."""
        return [
            SubRegion(FeatureLocation(prediction.start, prediction.end),
                      tool="boundary_fungal", label=prediction.anchor_id)
            for prediction in self.predictions
        ]


def count_motifs(promoter: str, linear_promoter: str,
                 capped_max_count: int | None = config.MAX_MOTIFS_PER_PROMOTER
                 ) -> utils.MotifCounts:
    """Counts the repeat and linear motifs in a gene's two promoters (including
        the linear promoter).

        Repeat motifs are counted over the full promoter and linear motifs
        over the shorter linear promoter, then the two sets are concatenated.
        The returned motif IDs are in increasing order.

        Arguments:
            promoter: the full promoter sequence, used for repeat motifs
            linear_promoter: the shorter promoter used for linear motifs
            capped_max_count: the maximum count kept for each motif in
                              either promoter, or None to keep uncapped counts

        Returns:
            a MotifCounts instance covering both kinds of motif
    """
    repeat_counts = utils.count_repeat_motifs(
        promoter, capped_max_count=capped_max_count)
    linear_counts = utils.count_linear_motifs(
        linear_promoter, capped_max_count=capped_max_count)

    # Repeat IDs precede linear IDs, preserving sorted order here.
    return utils.MotifCounts(
        motif_ids=np.concatenate((repeat_counts.motif_ids,
                                  linear_counts.motif_ids)),
        counts=np.concatenate((repeat_counts.counts, linear_counts.counts)),
    )


def get_motif_id_to_bg_freq(
        genes: dict[str, GeneWithMotifs], *,
        motif_type_id_boundaries: tuple[int, int] = config.motif_type_id_boundaries,
) -> dict[int, float]:
    """ motif ID mapped to its frequency within its own motif type

        Counts are summed over every gene, using the capped counts, and
        each motif's total is then divided by the total for its motif type.
        The frequencies of each observed type therefore sum to 1, so a motif
        can only be enriched relative to other motifs of the same kind.

        Arguments:
            genes: the contig's genes in CDS order, keyed by gene name
            motif_type_id_boundaries: the first direct-repeat and linear motif IDs,
                                      in that order

        Returns:
            each observed motif ID mapped to its background frequency
    """
    background_freqs: dict[int, float] = {}
    occurrences_per_motif_type = {
        motif_type: 0.0 for motif_type in utils.MotifType
    }

    # Sum counts for each motif and each motif type across all genes.
    for gene in genes.values():
        for motif_id, count in zip(gene.motif_ids, gene.capped_motif_counts):
            motif_id = int(motif_id)
            count = float(count)
            motif_type = utils.get_motif_type(
                motif_id, motif_type_id_boundaries=motif_type_id_boundaries)
            background_freqs[motif_id] = (
                background_freqs.get(motif_id, 0.0) + count)
            occurrences_per_motif_type[motif_type] += count

    # Divide each motif's count by the total count for its type.
    for motif_id in background_freqs:
        motif_type = utils.get_motif_type(
            motif_id, motif_type_id_boundaries=motif_type_id_boundaries)
        background_freqs[motif_id] /= occurrences_per_motif_type[motif_type]

    return background_freqs


@dataclass
class _WindowCandidates:
    """ Motif evidence for candidate windows sharing a starting gene.

        Gene indices refer to the anchor's neighbourhood
    """

    window_start_gene_index: int
    window_stop_gene_indices: np.ndarray
    observed: np.ndarray
    expected: np.ndarray
    eligible: np.ndarray
    in_anchor: np.ndarray
    in_adjacent: np.ndarray

    def select_best_motif_candidates(self, motif_ids: np.ndarray, anchor_id: str,
                                     anchor_index: int
                                     ) -> Iterator[EnrichedMotif]:
        """ Yields every eligible motif and window tied for the highest LLR.

            Scores eligible candidates using the log-likelihood ratio
            of observed to expected motif counts.

            Arguments:
                motif_ids: motif IDs
                anchor_id: the name of the anchor gene
                anchor_index: the anchor's index within its gene neighbourhood

            Returns:
                an iterator of EnrichedMotif instances tied for the highest LLR
                among these candidates, empty if no candidates are eligible
        """
        rows, columns = np.nonzero(self.eligible)
        if columns.size == 0:
            return
        hits = self.observed[rows, columns]
        expected_hits = self.expected[rows, columns]
        scores = hits * np.log(hits / expected_hits) - (hits - expected_hits)
        highest_llr = float(scores.max())
        for winner in np.flatnonzero(scores == highest_llr):
            row, column = rows[winner], columns[winner]
            if self.in_anchor[column]:
                mode = SelectionMode.ANCHOR
            elif self.in_adjacent[row, column]:
                mode = SelectionMode.ADJACENT
            else:
                mode = SelectionMode.ALL
            yield EnrichedMotif(
                anchor_id=anchor_id,
                motif_id=int(motif_ids[column]),
                llr=highest_llr,
                mode=mode,
                first_gene_offset=self.window_start_gene_index - anchor_index,
                gene_count=int(self.window_stop_gene_indices[row]) - self.window_start_gene_index,
            )


@dataclass
class _WindowMotifCounts:
    """ Cumulative motif counts for the genes surrounding one anchor.

        Each cumulative array starts with a zero row, so subtracting the
        row at a window's start from the row at its stop gives its total.
    """

    motif_ids: np.ndarray
    motif_types: np.ndarray
    background_frequencies: np.ndarray
    occurrences: np.ndarray
    genes_with_motif: np.ndarray
    type_occurrences: np.ndarray
    anchor_index: int

    @classmethod
    def calculate_cumulative_motif_counts(  # pylint: disable=too-many-locals
            cls, genes: list[GeneWithMotifs], anchor_index: int,
            motif_ids: np.ndarray, background_frequencies: dict[int, float], *,
            motif_type_id_boundaries: tuple[int, int] = config.motif_type_id_boundaries,
    ) -> _WindowMotifCounts:
        """ Calculates cumulative motif counts across genes surrounding an anchor.

            Computes running totals (cumulative) across genes so that any candidate window's
            counts can be calculated by subtracting its start from its stop.

            Arguments:
                genes: genes in the anchor's neighbourhood
                anchor_index: the position of the anchor within genes
                motif_ids: complete ordered range of motif IDs
                background_frequencies: each motif ID mapped to its contig-wide
                                        frequency within its motif type
                motif_type_id_boundaries: the first direct-repeat and linear
                                          motif IDs, in that order

            Returns:
                a _WindowMotifCounts instance with cumulative evidence for
                candidate window scoring
        """
        num_motif_types = len(utils.MotifType)
        # initialize windows length x num_motifs
        # which stores the count data, for the genes in the current window
        # that is gene_0 has 2 occurrences of motif_id 10
        counts = np.zeros((len(genes), len(motif_ids)), dtype=np.int64)
        # Store counts for each gene and motif type.
        type_counts = np.zeros((len(genes), num_motif_types), dtype=np.int64)
        # Loop through each gene object and its zero-based position within
        # the window.
        for index, gene in enumerate(genes):
            counts[index, gene.motif_ids] = gene.capped_motif_counts
            # To find out what kind of motif_type (inverted, direct, or linear) the motif_ids is
            gene_types = np.searchsorted(motif_type_id_boundaries,
                                         gene.motif_ids, side="right")
            # adds up motif-id occurrences by motif-type for the current gene
            # np.bincount: Count number of occurrences of each value in array
            type_counts[index] = np.bincount(
                gene_types, weights=gene.capped_motif_counts, minlength=num_motif_types)

        # initialize a 2d matrix: (number of genes in the window +1) x motif IDs
        # +1 extra row is used for cumulative sums later
        occurrences = np.zeros((len(genes) + 1, len(motif_ids)), dtype=np.int64)
        genes_with_motif = np.zeros_like(occurrences)
        # Store cumulative motif occurrence counts by type:
        # inverted, direct, linear.
        type_occurrences = np.zeros((len(genes) + 1, num_motif_types), dtype=np.int64)

        # Calculate running totals for each motif across genes.
        np.cumsum(counts, axis=0, out=occurrences[1:])
        np.cumsum(counts > 0, axis=0, out=genes_with_motif[1:])
        np.cumsum(type_counts, axis=0, out=type_occurrences[1:])
        return cls(
            motif_ids=motif_ids,
            # To find out what kind of motif_type (inverted, direct, or linear)
            # the motif_ids is
            motif_types=np.searchsorted(motif_type_id_boundaries, motif_ids, side="right"),
            background_frequencies=np.array([
                background_frequencies.get(motif_id, 0.0) for motif_id in motif_ids
            ], dtype=float),
            occurrences=occurrences,
            genes_with_motif=genes_with_motif,
            type_occurrences=type_occurrences,
            anchor_index=anchor_index,
        )

    def calculate_candidate_window_motif_counts(
            self, window_start_gene_index: int, min_genes_with_motif: int,
            max_genes_in_bgc: int) -> _WindowCandidates:
        """ Builds motif counts and eligibility masks for windows with the given start.

            Arguments:
                window_start_gene_index: starting gene index within the
                                        neighbourhood supplied to
                                        calculate_cumulative_motif_counts
                min_genes_with_motif: fewest genes whose promoters must contain
                                     the motif for it to be eligible
                max_genes_in_bgc: most genes a candidate window may span
        """
        first_stop = max(self.anchor_index + 1, window_start_gene_index + min_genes_with_motif)
        last_stop = min(len(self.occurrences) - 1, window_start_gene_index + max_genes_in_bgc)
        # List every allowed stop index.
        window_stops = np.arange(first_stop, last_stop + 1)

        # Subtract the cumulative counts before the start from the counts at each stop.
        observed = self.occurrences[window_stops] - self.occurrences[window_start_gene_index]
        genes_with_motif = (
            self.genes_with_motif[window_stops] - self.genes_with_motif[window_start_gene_index]
        )
        type_totals = (
            self.type_occurrences[window_stops] - self.type_occurrences[window_start_gene_index]
        )

        # The baseline uses each motif's contig-wide frequency within its type.
        expected = type_totals[:, self.motif_types] * self.background_frequencies
        enough_genes = genes_with_motif >= min_genes_with_motif
        enriched = observed > expected
        in_anchor, in_adjacent = self._motif_presence_near_anchor(
            window_start_gene_index, window_stops,
        )
        return _WindowCandidates(
            window_start_gene_index=window_start_gene_index,
            window_stop_gene_indices=window_stops,
            observed=observed,
            expected=expected,
            eligible=enough_genes & enriched,
            in_anchor=in_anchor,
            in_adjacent=in_adjacent,
        )

    def _motif_presence_near_anchor(
            self, window_start: int, window_stops: np.ndarray,
    ) -> tuple[np.ndarray, np.ndarray]:
        """ Finds motifs in the anchor and its immediate neighbourhood.

            The neighbourhood includes the anchor and is clipped to each candidate window.
        """
        anchor = self.anchor_index
        presence_counts = self.genes_with_motif
        in_anchor = (presence_counts[anchor + 1] - presence_counts[anchor]) > 0

        nearby_start = max(window_start, anchor - 1)
        nearby_stops = np.minimum(window_stops, anchor + 2)
        in_adjacent = (presence_counts[nearby_stops] - presence_counts[nearby_start]) > 0
        return in_anchor, in_adjacent


def _rank_motif_candidate(result: EnrichedMotif, *,
                          linear_motif_id_start: int = config.linear_motif_id_start
                          ) -> tuple[float, bool, SelectionMode, int, int]:
    """ Builds a comparison key for a motif candidate.

        Ranks by LLR, repeat preference, proximity, window size, and motif ID.

        Arguments:
            result: the motif candidate to rank
            linear_motif_id_start: the first ID reserved for linear motifs

        Returns:
            a comparison key where smaller values indicate better candidates
    """
    # Break equal-LLR ties using the remaining ranking criteria.
    return (
        -result.llr,                                    # Higher LLR first.
        result.motif_id >= linear_motif_id_start,       # Repeats win when LLR scores tie.
        result.mode,                                    # Anchor, then adjacent, then all.
        result.gene_count,                              # Fewer genes in the window first.
        result.motif_id,                                # Lower motif ID breaks the final tie.
    )


def _score_anchor(neighbourhood_counts: _WindowMotifCounts, anchor_id: str,
                  min_genes_with_motif: int, max_genes_in_bgc: int, *,
                  linear_motif_id_start: int = config.linear_motif_id_start
                  ) -> EnrichedMotif | None:
    """ Selects the best eligible motif window containing the anchor.

        Arguments:
            neighbourhood_counts: cumulative motif counts for the anchor
                                  and surrounding genes
            anchor_id: the name of the anchor gene
            min_genes_with_motif: minimum number of genes whose promoters
                                  must contain the motif within a window
            max_genes_in_bgc: maximum number of genes in a candidate window
            linear_motif_id_start: first linear motif ID, used to distinguish
                                   linear motifs from repeat motifs

        Returns:
            the best candidate, or None if no eligible motif window exists
    """
    best = None
    best_rank = None
    # Try different starting genes for candidate windows containing the anchor gene.
    # The maximum window size limits how far before the anchor a window can start.
    for window_start_gene_index in range(
            max(0, neighbourhood_counts.anchor_index - max_genes_in_bgc + 1),
            neighbourhood_counts.anchor_index + 1):
        candidates = neighbourhood_counts.calculate_candidate_window_motif_counts(
            window_start_gene_index, min_genes_with_motif, max_genes_in_bgc)
        for candidate in candidates.select_best_motif_candidates(
                neighbourhood_counts.motif_ids, anchor_id,
                neighbourhood_counts.anchor_index):
            candidate_rank = _rank_motif_candidate(
                candidate, linear_motif_id_start=linear_motif_id_start)
            # Keep the first candidate when all ranking values are equal.
            if best_rank is None or candidate_rank < best_rank:
                best_rank = candidate_rank
                best = candidate
    return best


def find_enriched_motifs(  # pylint: disable=too-many-locals
        genes: dict[str, GeneWithMotifs], background_frequencies: dict[int, float],
        min_genes_with_motif: int = config.MIN_GENES_WITH_MOTIF,
        max_genes_in_bgc: int = config.MAX_GENES_IN_BGC,
        *, num_motifs: int = config.n_motifs,
        motif_type_id_boundaries: tuple[int, int] = config.motif_type_id_boundaries,
) -> dict[str, EnrichedMotif | None]:
    """ Finds the most enriched motif in a window of genes around each anchor

        Windows are ranked by a log-likelihood ratio (LLR) of observed
        against expected counts, the expectation being the motif's background
        frequency times the total count of its motif type in the window. To be
        eligible a motif must occur in at least min_genes_with_motif genes of
        the window and exceed its expected count.

        Arguments:
            genes: the contig's genes in CDS order, keyed by gene name
            background_frequencies: each motif ID mapped to its background
                                    frequency within its own motif type
            min_genes_with_motif: the fewest genes of a window a motif must
                                  occur in to be considered
            max_genes_in_bgc: the most genes a window may span, and the
                              furthest it may reach from its anchor
            num_motifs: the number of reserved motif IDs
            motif_type_id_boundaries: the first direct-repeat and linear motif IDs,
                                      in that order

        Returns:
            each core anchor's name mapped to its best EnrichedMotif, or to
            None where no motif was eligible
    """
    motif_ids = np.arange(num_motifs)
    gene_list = list(genes.values())
    results: dict[str, EnrichedMotif | None] = {}
    for anchor_index, anchor in enumerate(gene_list):
        if not anchor.is_anchor:
            continue
        # Keep the anchor neighbourhood within the contig.
        first_index = max(0, anchor_index - max_genes_in_bgc)
        stop_index = min(len(gene_list), anchor_index + max_genes_in_bgc + 1)
        neighbourhood_counts = _WindowMotifCounts.calculate_cumulative_motif_counts(
            gene_list[first_index:stop_index], anchor_index - first_index,
            motif_ids, background_frequencies,
            motif_type_id_boundaries=motif_type_id_boundaries,
        )
        results[anchor.name] = _score_anchor(
            neighbourhood_counts, anchor.name, min_genes_with_motif,
            max_genes_in_bgc,
            linear_motif_id_start=motif_type_id_boundaries[1],
        )
    return results


def build_report(contig: str, genes: dict[str, GeneWithMotifs],
                 enriched_motifs: dict[str, EnrichedMotif | None], *,
                 min_llr: float = config.MIN_LLR
                 ) -> BoundaryDetectionResults:
    """ Collects the accepted results into a report for one contig.

        Only results with an LLR of at least min_llr are
        reported. Each accepted window is located from its anchor's
        position plus the result's first_gene_offset, spanning gene_count
        genes.

        Arguments:
            contig: the record ID the report covers
            genes: the contig's genes in CDS order, keyed by gene name
            enriched_motifs: each anchor's result
            min_llr: the inclusive minimum LLR for reported boundaries

        Returns:
            the accepted boundary predictions
    """
    gene_list = list(genes.values())
    gene_indices = {gene.name: index for index, gene in enumerate(gene_list)}
    predictions: list[BoundaryPrediction] = []
    for motif in enriched_motifs.values():
        if motif is None or motif.llr < min_llr:
            continue

        first_index = gene_indices[motif.anchor_id] + motif.first_gene_offset
        window_genes = gene_list[first_index:first_index + motif.gene_count]
        motif_type = utils.get_motif_type(motif.motif_id)
        motif_counts = {
            gene.name: int(gene.uncapped_motif_counts[
                gene.motif_ids == motif.motif_id].sum())
            for gene in window_genes
        }

        predictions.append(BoundaryPrediction(
            anchor_id=motif.anchor_id,
            start=min(gene.start for gene in window_genes),
            end=max(gene.end for gene in window_genes),
            motif_type=motif_type.value,
            motif_string_dna=utils.convert_motif_id_to_text(motif.motif_id),
            motif_counts=motif_counts,
        ))
    return BoundaryDetectionResults(contig=contig, predictions=predictions)


def predict_boundaries(record: Record) -> BoundaryDetectionResults:
    """ Predict fungal BGC boundaries from promoter motifs around core CDS features.

        Compare observed motif counts in candidate windows with counts expected
        from each motif's contig-wide frequency within its motif type.

        Arguments:
            record: the record to analyse

        Returns:
            a BoundaryDetectionResults object. Its predictions list is empty
            when the contig has fewer than config.MIN_GENES_PER_CONTIG CDS
            features or no eligible window meets the minimum log-likelihood ratio.
    """

    cds_features = record.get_cds_features()
    if len(cds_features) < config.MIN_GENES_PER_CONTIG:
        return BoundaryDetectionResults(contig=record.id, predictions=[])

    annotated_genes = {}
    for cds in cds_features:
        gene = GeneWithMotifs.from_feature(record, cds)
        annotated_genes[gene.name] = gene

    motif_id_to_bg_freq = get_motif_id_to_bg_freq(annotated_genes)

    enriched_motifs = find_enriched_motifs(
        genes=annotated_genes,
        background_frequencies=motif_id_to_bg_freq,
        max_genes_in_bgc=config.MAX_GENES_IN_BGC,
    )

    return build_report(record.id, annotated_genes, enriched_motifs)
