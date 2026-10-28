# License: GNU Affero General Public License v3 or later
# A copy of GNU AGPL v3 should have been included in this software package in
# LICENSE.txt.
"""Predict fungal cluster boundaries from motif enrichment in promoters."""


from dataclasses import dataclass
from enum import StrEnum, auto
import random
import re

from Bio.Seq import reverse_complement
import numpy as np

from antismash.common.secmet import CDSFeature, Record
from antismash.modules.clusterblast.core import get_core_gene_ids

from . import utils
from . import config


class SelectionMode(StrEnum):
    """Which promoters may supply motifs for an anchor."""

    ANCHOR = auto()
    ADJACENT = auto()
    ALL = auto()


class MotifType(StrEnum):
    """The DNA arrangement represented by a motif."""

    INVERTED = auto()
    DIRECT = auto()
    LINEAR = auto()


# Data Structures
@dataclass
class GeneWithPromoters:  # pylint: disable=too-many-instance-attributes
    """One CDS feature's annotation and promoter sequences for motif searches.

    Both promoters follow the gene's reading direction, with position 0 at the
    upstream end of the extracted sequence. Each distance counts bases from
    that end to the annotated CDS start within its corresponding promoter;
    this does not identify the actual transcription start site.
    """
    name: str
    start: int
    end: int
    strand: int
    promoter: str
    promoter_to_cds_start_dist: int
    linear_promoter: str
    linear_promoter_to_cds_start_dist: int

    @classmethod
    def from_feature(cls, record: Record,
                     cds: CDSFeature) -> "GeneWithPromoters":
        """Extract promoter windows for a CDS.

            Include up to config.UPSTREAM upstream bases and config.DOWNSTREAM
            downstream bases.
            Reverse-complement minus strand windows, so all promoters has the same direction.
            Linear window (used for only linear motif searches later) has its own
            config.LINEAR_UPSTREAM upstream bases and the same downstream sequence.

            Arguments:
                record: the record supplying the sequence containing the CDS
                cds: the CDS feature supplying the gene name, location
                and strand

            Returns:
                A new instance with the gene name, zero-based start and
                end-exclusive genomic coordinates, strand, and the full and
                linear promoter sequences. The promoter_to_cds_start_dist and
                linear_promoter_to_cds_start_dist fields measure the distance
                in bases from the start of their respective promoter sequences
                to the annotated 5' CDS boundary.
        """
        cds_start_coord = int(cds.location.start)
        cds_end_coord = int(cds.location.end)
        strand = cds.location.strand

        # Keep the promoter window within the contig.
        if strand == -1:
            promoter_window_start = max(0, cds_end_coord - config.DOWNSTREAM)
            promoter_window_end = min(len(record.seq), cds_end_coord + config.UPSTREAM)
            promoter_seq = str(
                record.seq[promoter_window_start:promoter_window_end].reverse_complement()
            )
            promoter_to_cds_start_dist = promoter_window_end - cds_end_coord
        else:
            promoter_window_start = max(0, cds_start_coord - config.UPSTREAM)
            promoter_window_end = min(len(record.seq), cds_start_coord + config.DOWNSTREAM)
            promoter_seq = str(record.seq[promoter_window_start:promoter_window_end])
            promoter_to_cds_start_dist = cds_start_coord - promoter_window_start

        # Linear motifs are searched closer to the gene than repeats are.
        crop = max(0, promoter_to_cds_start_dist - config.LINEAR_UPSTREAM)
        return cls(name=cds.get_name(),
                   start=cds_start_coord, end=cds_end_coord, strand=strand,
                   promoter=promoter_seq, promoter_to_cds_start_dist=promoter_to_cds_start_dist,
                   linear_promoter=promoter_seq[crop:],
                   linear_promoter_to_cds_start_dist=promoter_to_cds_start_dist - crop)


@dataclass
class Contig:
    """One contig's genes, promoter motif counts and core gene anchors.

    genes and gene_motif_counts are parallel lists: entry i of each
    describes the same gene. The integers in anchors index into genes.
    motif_totals describes the total capped number of each motif across the
    contig from gene_motif_counts.
    """
    id: str
    genes: list[GeneWithPromoters]
    gene_motif_counts: list[utils.MotifCounts]
    motif_totals: np.ndarray
    anchors: list[tuple[str, int]]


@dataclass
class AnchorScores:  # pylint: disable=too-many-instance-attributes
    """The best window found for every motif around one anchor.

    Every array but neighbourhood_counts is indexed by motif ID, and
    describes that motif's highest-scoring window. The gene offsets are
    relative to the anchor; the other fields locate the anchor itself.
    """
    anchor_index: int              # the anchor's index in the contig's genes
    center: int                    # the anchor's row in neighbourhood_counts
    neighbourhood_counts: np.ndarray  # motif counts per neighbouring gene
    llr: np.ndarray                # enrichment score of the best window
    observed_hits: np.ndarray      # hits in that window (capped counts)
    genes_with_motif: np.ndarray   # genes in that window carrying the motif
    expected_hits: np.ndarray      # hits the background rate predicts there
    first_gene_offset: np.ndarray  # window start, in genes from the anchor
    last_gene_offset: np.ndarray   # window end, in genes from the anchor
    gene_counts: np.ndarray        # genes in that window


@dataclass
class GeneRow:  # pylint: disable=too-many-instance-attributes
    """One gene in the evidence table for a chosen motif window.

    The table covers the window plus config.BUFFER_GENES genes on each
    side, clipped at the contig ends, so some rows lie outside the
    prediction; in_cluster marks the ones inside it.

    start, end and position are in contig order and say nothing about the
    gene's strand. hit_positions_from_gene_start is in promoter
    coordinates, which follow the gene's own 5' to 3' direction, so
    negative means upstream for forward and reverse genes alike.
    """
    name: str
    position: int                              # genes from the anchor
    in_cluster: bool                           # inside the chosen window
    start: int
    end: int
    strand: int
    motif_hits: int                            # uncapped number of occurrences
    hit_positions_from_gene_start: list[int]   # negative = upstream


@dataclass
class Candidate:  # pylint: disable=too-many-instance-attributes
    """A motif enriched around an anchor, together with its best window."""
    mode: SelectionMode
    motif_type: MotifType
    motif: str                  # readable form, e.g. "ACG-N2-CGT"
    motif_id: int
    first_gene_offset: int      # inclusive window start, in genes from the anchor
    last_gene_offset: int       # inclusive window end, in genes from the anchor
    gene_count: int             # genes in the window
    llr: float                  # enrichment score
    llr_per_gene: float
    observed_hits: int          # hits in the window (capped counts)
    genes_with_motif: int       # genes in the window with the motif
    expected_sites: float       # hits the background predicts
    p_empirical: float          # how often random genes scored as high


@dataclass
class BoundaryPrediction:
    """An accepted boundary and the promoter evidence used for reporting.

    Gene IDs cover the chosen window. The gene table and coordinates also
    include the configured buffer genes around that window.
    """

    anchor_id: str
    motif_string_dna: str
    motif_type: MotifType
    p_empirical: float
    ids_of_genes_within_the_predicted_boundary: list[str]
    gene_table: list[GeneRow]

    @property
    def start(self) -> int:
        """ Return the boundary's inclusive start, including buffer genes.

            Arguments:
                none

            Returns:
                The smallest zero-based start coordinate among all genes in
                gene_table, relative to the contig sequence.
        """
        return min(gene.start for gene in self.gene_table)

    @property
    def end(self) -> int:
        """ Return the boundary's exclusive end, including buffer genes.

            Arguments:
                none

            Returns:
                The largest zero-based, end-exclusive coordinate among all
                genes in gene_table, relative to the contig sequence.
        """
        return max(gene.end for gene in self.gene_table)


@dataclass
class FundariesReport:
    """All accepted boundaries for one contig."""

    contig: str
    predictions: list[BoundaryPrediction]

# Read data and motif counting


def count_motifs(gene: GeneWithPromoters) -> utils.MotifCounts:
    """Return sorted motif IDs and capped counts for this gene's promoters.

    Arguments:
        gene: the gene supplying the repeat and linear promoter sequences

    Returns:
        Sorted unique motif IDs and counts capped at config.MAX_COUNTS.
    """

    repeat_counts = utils.repeat_motif_counts(gene.promoter)
    linear_counts = utils.linear_motif_counts(gene.linear_promoter)
    # Repeat IDs precede linear IDs, so concatenation preserves sorted order.
    return utils.MotifCounts(
        motif_ids=np.concatenate((repeat_counts.motif_ids, linear_counts.motif_ids)),
        counts=np.concatenate((repeat_counts.counts, linear_counts.counts)),
    )


def sum_motif_counts(gene_motif_counts: list[utils.MotifCounts]) -> np.ndarray:
    """Sum each motif's capped counts across the contig.

    Arguments:
        gene_motif_counts: each gene's motif IDs and their capped counts

    Returns:
        One entry per motif ID, holding that motif's capped count summed over
        every gene on the contig.
    """
    motif_totals = np.zeros(config.N_MOTIFS)
    for gene_counts in gene_motif_counts:
        motif_totals[gene_counts.motif_ids] += gene_counts.counts
    return motif_totals


# Scoring functions

def get_neighbourhood_indices(genes: list[GeneWithPromoters],
                              index: int, reach: int) -> range:
    """Returns neighbourhood of genes from the given anchor"""
    first_index = max(0, index - reach)
    stop_index = min(len(genes), index + reach + 1)
    return range(first_index, stop_index)


def enumerate_gene_windows(n_genes: int,
                           center: int) -> tuple[np.ndarray, np.ndarray]:
    """Returns every window of neighbouring genes that contains the anchor."""
    window_starts = []
    window_ends = []
    for start in range(max(0, center - config.MAX_GENES + 1), center + 1):
        for end in range(center, min(n_genes, start + config.MAX_GENES)):
            if end - start + 1 >= config.MIN_GENES:
                window_starts.append(start)
                window_ends.append(end)
    return np.array(window_starts), np.array(window_ends)


def score_windows(genes: list[GeneWithPromoters],
                  gene_motif_counts: list[utils.MotifCounts],
                  motif_totals: np.ndarray, anchor_index: int) -> AnchorScores:
    """For every motif, find the window of genes around this one where it is
    most enriched compared to the rest of the contig.

    A window scores high when the motif occurs there more often than the
    background rate predicts. The score is a Poisson log-likelihood ratio:
    observed * log(observed / expected) - (observed - expected).

    Only the config.N_FLANKING genes either side of the anchor are examined,
    clipped at the contig ends, and only windows that contain the anchor and
    hold between config.MIN_GENES and config.MAX_GENES genes are scored. The
    background rate is each motif's count per gene across the genes outside
    that neighbourhood, with a floor of half a hit among them, so a motif
    never seen outside still has something to be compared against.

    Every motif gets an entry, including motifs that are never enriched:
    those keep a score of 0 and the first window enumerated, so check llr
    before reading a motif's window.

    Arguments:
        genes: all the genes in the contig, used to locate the neighbourhood
            and to count the genes outside it
        gene_motif_counts: each gene's motif IDs and their capped counts, in
            the same order as genes
        motif_totals: total count of each motif across the contig, from which
            the background rate is taken
        anchor_index: the anchor's index in genes

    Returns:
        The AnchorScores for this anchor: per-motif arrays describing each
        motif's best window, plus the anchor's index, its row in the count
        table (center) and the table itself (neighbourhood_counts).
    """
    # Count every motif in each gene near the anchor
    neighbour_indices = get_neighbourhood_indices(genes, anchor_index,
                                                  config.N_FLANKING)
    # Center represent also the anchor gene index, inside of the window,
    # but recalculate if it is not in the center
    center = anchor_index - neighbour_indices[0]

    neighbourhood_counts = np.zeros((len(neighbour_indices), config.N_MOTIFS))
    for row, gene_counts in enumerate(
        gene_motif_counts[neighbour_indices.start:neighbour_indices.stop]
    ):
        neighbourhood_counts[row, gene_counts.motif_ids] = gene_counts.counts

    # Calculate the background rate
    genes_outside = len(genes) - len(neighbour_indices)
    if genes_outside < 1:
        raise ValueError(f"This contig leaves no promoters "
                         f"outside the {len(neighbour_indices)}-gene "
                         f"neighbourhood to compare against")

    rate = (motif_totals - neighbourhood_counts.sum(axis=0)) / genes_outside
    rate = np.maximum(rate, 0.5 / genes_outside)

    # Score possible winners together for each left endpoint.
    window_starts, window_ends = enumerate_gene_windows(len(neighbour_indices), center)
    best_windows = utils.find_best_windows(neighbourhood_counts, rate,
                                           window_starts, window_ends, center)
    best_window_size = best_windows.right - best_windows.left + 1

    return AnchorScores(
        anchor_index=anchor_index,
        center=center,
        neighbourhood_counts=neighbourhood_counts,
        llr=best_windows.llr,
        observed_hits=best_windows.observed_hits,
        genes_with_motif=best_windows.genes_with_motif,
        expected_hits=best_window_size * rate,
        first_gene_offset=best_windows.left - center,
        last_gene_offset=best_windows.right - center,
        gene_counts=best_window_size,
    )


def select_motifs_for_mode(anchor_scores: AnchorScores,
                           mode: SelectionMode) -> np.ndarray:
    """IDs of the motifs a mode may choose from.

    "anchor":   motifs present in the anchor gene's own promoter
    "adjacent": motifs present in the anchor or the genes either side of it
    "all":      every motif

    Arguments:
        anchor_scores: one anchor's scores, for the neighbourhood counts and
            anchor's row in them
        mode: a SelectionMode member

    Returns:
        Motif IDs of the different modes.

    Raises:
        ValueError: when mode is not a recognised selection mode.
    """
    neighbourhood_counts = anchor_scores.neighbourhood_counts
    center = anchor_scores.center

    if mode == SelectionMode.ANCHOR:
        return np.flatnonzero(neighbourhood_counts[center])
    if mode == SelectionMode.ADJACENT:
        adj_genes = neighbourhood_counts[max(0, center - 1):center + 2]
        return np.flatnonzero(adj_genes.sum(axis=0))
    if mode == SelectionMode.ALL:
        return np.arange(config.N_MOTIFS)
    raise ValueError(f"Unknown motif selection mode: {mode!r}")


def find_highest_llr(anchor_scores: AnchorScores) -> dict[tuple[SelectionMode, MotifType], float]:
    """Return the highest score reached in each (mode, motif type) category.

    A mode restricts which motifs may be chosen and a motif type is a block of
    the ID range, so the nine combinations are nine separate searches over the
    same scores.

    Arguments:
        anchor_scores: one anchor's per-motif scores, from score_windows

    Returns:
        One entry for each (mode, motif type) pair, holding the highest score
        among the motifs that the mode may choose from and that belong to that
        type. An entry is 0.0 when the mode offers no motif of the type, or
        when none of them scored.
    """
    motif_llr = anchor_scores.llr
    best_per_category = {}
    for mode in SelectionMode:
        mode_motif_ids = select_motifs_for_mode(anchor_scores, mode)
        for motif_type in MotifType:
            type_start, type_stop = config.MOTIF_TYPE_RANGE[motif_type]
            # Motifs of this mode that belong to this type.
            in_type = (mode_motif_ids >= type_start) & (mode_motif_ids < type_stop)
            type_motif_ids = mode_motif_ids[in_type]
            # Highest score among them, or 0 if there are none.
            if len(type_motif_ids) == 0:
                top_score = 0.0
            else:
                top_score = float(motif_llr[type_motif_ids].max())
            best_per_category[(mode, motif_type)] = top_score
    return best_per_category


# calibration against random genes

def calibrate(genes: list[GeneWithPromoters],
              gene_motif_counts: list[utils.MotifCounts],
              motif_totals: np.ndarray) -> dict[tuple[SelectionMode, MotifType], list[float]]:
    """Run the same window search at random genes to see what chance produces.

    Thousands of motif and window combinations are tried per anchor, so a high
    score means little on its own: it has to be compared against the high
    scores reached at genes picked at random from this contig, which carries
    its own GC content and gene density. Only the best score per (mode, motif
    type) is kept for each random gene, so the comparison already accounts for
    all those combinations having been tried.

    The sample is drawn from the genes at least config.MIN_FLANK from either
    end of the contig, so every random anchor has neighbours on both sides,
    and holds up to config.NUM_TEST_ANCHORS of them. config.SEED fixes the
    draw, so repeated runs on the same contig calibrate identically.

    Arguments:
        genes: all the genes contained in the contig
        gene_motif_counts: each gene's motif IDs and their counts, in the same
            order as genes
        motif_totals: total count of each motif across the contig, used for
            the background rate

    Returns:
        One entry for each (mode, motif type) pair, holding the best score
        reached at each random gene: one score per sampled gene, 0.0 where no
        motif of that type scored.

    Raises:
        ValueError: when the contig has too few genes to sample from.
    """
    selected_anchors_indices = range(config.MIN_FLANK, len(genes) - config.MIN_FLANK)
    if not selected_anchors_indices:
        raise ValueError("Too few genes to calibrate on")
    rng = random.Random(config.SEED)
    sample_size = min(config.NUM_TEST_ANCHORS, len(selected_anchors_indices))
    random_genes_indices = rng.sample(selected_anchors_indices, sample_size)

    calibration_scores: dict[tuple[SelectionMode, MotifType], list[float]] = {
        (mode, motif_type): []
        for mode in SelectionMode for motif_type in MotifType
    }

    for anchor_index in random_genes_indices:
        anchor_scores = score_windows(genes, gene_motif_counts, motif_totals,
                                      anchor_index)
        best_per_category = find_highest_llr(anchor_scores)

        for category, best_score in best_per_category.items():
            calibration_scores[category].append(best_score)

    return calibration_scores


def calculate_empirical_p_value(llr: float, calibration_scores: list[float]) -> float:
    """How often chance matched this score, counting the observation itself.

    Arguments:
        llr: the candidate's enrichment score
        calibration_scores: the scores random genes reached in this
            candidate's own (mode, motif type) category

    Returns:
        The fraction of those scores that matched or beat llr, from
        1 / (len(calibration_scores) + 1) up to 1.
    """
    times_higher = sum(value >= llr for value in calibration_scores)
    p_empirical = (1 + times_higher) / (len(calibration_scores) + 1)
    return p_empirical


# Decision module: build, check and rank results
def build_candidates(
    anchor_scores: AnchorScores,
    calibration_scores: dict[tuple[SelectionMode, MotifType], list[float]],
) -> list[Candidate]:
    """Build candidates from all the possible motifs. One candidate per mode
    and enriched motif.

    Only motifs that scored above zero are kept. The rest never beat their
    background, and they carry whichever window happened to be tried first,
    which means nothing.

    The modes are nested, so the same motif can appear up to three times, once
    per mode that offers it. Each copy gets its own empirical p-value from
    that mode's calibration, because a motif sitting in the anchor's own
    promoter is a rarer find than the same motif somewhere in the
    neighbourhood, and the two should not be judged on one scale.

    Counts are converted to plain ints and floats so that no numpy scalar
    reaches the report.

    Arguments:
        anchor_scores: one anchor's per-motif scores and windows, from
            score_windows
        calibration_scores: the scores random genes reached in each (mode,
            motif type) category, from calibrate

    Returns:
        One Candidate per mode and enriched motif, in mode order. Empty when
        no motif scored at this anchor.
    """
    candidates = []
    for mode in SelectionMode:
        motif_ids = select_motifs_for_mode(anchor_scores, mode)
        for motif_id in motif_ids[anchor_scores.llr[motif_ids] > 0]:
            motif_type = MotifType(utils.motif_type(motif_id))
            llr = float(anchor_scores.llr[motif_id])
            gene_count = int(anchor_scores.gene_counts[motif_id])
            candidates.append(Candidate(
                mode=mode,
                motif_type=motif_type,
                motif=utils.motif_int_to_text(motif_id),
                motif_id=int(motif_id),
                first_gene_offset=int(anchor_scores.first_gene_offset[motif_id]),
                last_gene_offset=int(anchor_scores.last_gene_offset[motif_id]),
                gene_count=gene_count,
                llr=llr,
                llr_per_gene=llr / gene_count,
                observed_hits=int(anchor_scores.observed_hits[motif_id]),
                genes_with_motif=int(anchor_scores.genes_with_motif[motif_id]),
                expected_sites=float(anchor_scores.expected_hits[motif_id]),
                p_empirical=calculate_empirical_p_value(llr, calibration_scores[mode, motif_type]),
            ))
    return candidates


def passes_checks(candidate: Candidate) -> bool:
    """Whether a candidate passes every acceptance rule.

    The rules are independent, so the first failure settles it: the candidate
    is not reported and nothing records which rule stopped it.

    Arguments:
        candidate: a boundary candidate

    Returns:
        True when the candidate clears every rule in config, False as soon as
        one fails.
    """
    if candidate.p_empirical >= config.P_VAL_MAX:
        return False
    if candidate.genes_with_motif < config.MIN_GENES_WITH_MOTIFS:
        return False
    if candidate.gene_count > config.MAX_PREDICTED:
        return False
    if candidate.llr_per_gene < config.MIN_LLR_PER_GENE:
        return False
    support = candidate.genes_with_motif / candidate.gene_count
    if candidate.motif_type == MotifType.LINEAR and support < config.MIN_LINEAR_SUPPORT:
        return False
    return True


def get_candidate_sort_key(candidate: Candidate) -> tuple[bool, float, float, int, MotifType, int]:
    """Sort key for choosing the best candidate: smaller is better.

    The fields are a list of preferences, not a score: min() compares them in
    order and the first one that differs decides, so a repeat motif beats any
    linear motif whatever its p-value. The score is megative because it is the
    one field where larger is better. The last two options are just made for
    untie same candidates.

    Arguments:
        candidate: a candidate that passed every acceptance check

    Returns:
        A tuple to order candidates by, smallest first.
    """
    is_linear = candidate.motif_type == MotifType.LINEAR
    p_empirical = candidate.p_empirical
    negative_score = -candidate.llr
    mode_order = tuple(SelectionMode).index(candidate.mode)
    motif_type = candidate.motif_type
    motif_id = candidate.motif_id

    return (is_linear, p_empirical, negative_score, mode_order,
            motif_type, motif_id)


# Reporting functions

def build_motif_patterns(candidate: Candidate) -> set[str]:
    """Every spelling of a candidate's motif, as regular expressions.

    The spacer becomes dots, so "ACG-N2-CGT" turns into "ACG..CGT". A linear
    motif and a direct repeat each need their reverse complement as well; an
    inverted repeat does not, because it reads the same from either strand.

    Arguments:
        candidate: the candidate whose motif is to be matched

    Returns:
        Patterns(depending on motif type), for re.finditer to search a promoter with.
    """
    if candidate.motif_type == MotifType.LINEAR:
        word = candidate.motif
        return {word, reverse_complement(word)}
    arm, spacer, other = candidate.motif.split("-")
    gap = "." * int(spacer[1:])
    patterns = {arm + gap + other}
    if candidate.motif_type == MotifType.DIRECT:
        patterns.add(reverse_complement(arm) + gap +
                     reverse_complement(other))
    return patterns


def find_motif_positions(sequence: str, patterns: set[str],
                         tss_offset: int) -> list[int]:
    """Returns positions of every occurrence of the motif, relative
    to the gene start.

    Arguments:
        sequence: the promoter to search
        patterns: the motif's patterns, from build_motif_patterns
        tss_offset: where the annotated CDS start sits in sequence

    Returns:
        Offsets from the annotated CDS start, ascending, negative upstream.
    """
    hit_positions = []
    for pattern in patterns:
        matches = re.finditer(f"(?={pattern})", sequence)
        hit_positions += [match.start() - tss_offset for match in matches]
    return sorted(hit_positions)


def build_gene_table(genes: list[GeneWithPromoters], anchor_scores: AnchorScores,
                     candidate: Candidate) -> list[GeneRow]:
    """Returns a summary of where the motif occurs in each gene in its
    window.

    The table also covers config.BUFFER_GENES genes either side of the window,
    clipped at the contig ends, so the report shows what lies just outside the
    boundary. Linear motifs are searched in the shorter linear promoter and
    repeats in the full one, and these hits are uncapped, unlike the counts
    used for scoring.

    Arguments:
        genes: all of the contig's genes
        anchor_scores: the anchor's scores, for its index in genes
        candidate: the accepted candidate whose window is described

    Returns:
        One GeneRow per gene in the window and its margin, in contig order.
    """
    patterns = build_motif_patterns(candidate)
    anchor_index = anchor_scores.anchor_index
    # The margin stops at the end of the contig.
    first_index = max(0, anchor_index + candidate.first_gene_offset
                      - config.BUFFER_GENES)
    last_index = min(len(genes) - 1,
                     anchor_index + candidate.last_gene_offset
                     + config.BUFFER_GENES)
    rows = []
    for index in range(first_index, last_index + 1):
        gene = genes[index]
        if candidate.motif_type == MotifType.LINEAR:
            hit_positions = find_motif_positions(
                gene.linear_promoter, patterns,
                gene.linear_promoter_to_cds_start_dist)
        else:
            hit_positions = find_motif_positions(
                gene.promoter, patterns, gene.promoter_to_cds_start_dist)

        position = index - anchor_index
        rows.append(GeneRow(
            name=gene.name,
            position=position,
            in_cluster=candidate.first_gene_offset <= position <= candidate.last_gene_offset,
            start=gene.start,
            end=gene.end,
            strand=gene.strand,
            motif_hits=len(hit_positions),
            hit_positions_from_gene_start=hit_positions,
        ))
    return rows


def decide_boundaries(
    contig: Contig,
    calibration_scores: dict[tuple[SelectionMode, MotifType], list[float]],
) -> FundariesReport:
    """ Predict cluster boundaries around the contig's core genes.

        Score motif windows for each core gene and keep the best candidate
        that passes all acceptance checks.

        Arguments:
            contig: the contig supplying genes, promoter motif counts and core gene anchors
            calibration_scores: best scores from random genes for each
                                (search mode (anchor, adjacent, all), motif type),
                                used to compute empirical p-values

        Returns:
            A report containing the contig ID and accepted boundary predictions.
            Each prediction includes the selected motif, its empirical p-value,
            the genes within the boundary and a table showing motif hits for
            genes within and around the predicted boundary.
    """

    report = FundariesReport(contig=contig.id, predictions=[])
    for core_gene, core_gene_index in contig.anchors:

        anchor_scores = score_windows(contig.genes, contig.gene_motif_counts,
                                      contig.motif_totals, core_gene_index)
        candidates = build_candidates(anchor_scores, calibration_scores)

        passing = [candidate for candidate in candidates if passes_checks(candidate)]
        if not passing:
            continue

        best_candidate = min(passing, key=get_candidate_sort_key)
        first_gene_index = core_gene_index + best_candidate.first_gene_offset
        last_gene_index = core_gene_index + best_candidate.last_gene_offset
        report.predictions.append(BoundaryPrediction(
            anchor_id=core_gene,
            motif_string_dna=best_candidate.motif,
            motif_type=best_candidate.motif_type,
            p_empirical=best_candidate.p_empirical,
            ids_of_genes_within_the_predicted_boundary=[
                contig.genes[i].name for i in range(first_gene_index, last_gene_index + 1)
            ],
            gene_table=build_gene_table(contig.genes, anchor_scores, best_candidate),
        ))
    return report


def main(record: Record) -> FundariesReport:
    """Run Fundaries and report boundary predictions for a contig.

    Arguments:
        record: antiSMASH Record representing the contig

    Returns:
        A report containing the contig ID and its accepted boundary
        predictions.
    """
    # Step 1: read contig, find and count motifs.
    genes_with_promoters = [
        GeneWithPromoters.from_feature(record, cds)
        for cds in record.get_cds_features()
    ]
    gene_motif_counts = [count_motifs(gene) for gene in genes_with_promoters]
    motif_totals = sum_motif_counts(gene_motif_counts)
    core_gene_ids = get_core_gene_ids(record)
    anchors = []
    for index, gene in enumerate(genes_with_promoters):
        if gene.name in core_gene_ids:
            anchors.append((gene.name, index))

    contig = Contig(id=record.id, genes=genes_with_promoters,
                    gene_motif_counts=gene_motif_counts,
                    motif_totals=motif_totals, anchors=anchors)

    # Step 2: Find the best scores random genes reach by chance
    calibration_scores = calibrate(contig.genes, contig.gene_motif_counts, contig.motif_totals)

    # Step 3: Collect accepted boundaries for all core genes.
    return decide_boundaries(contig, calibration_scores)
