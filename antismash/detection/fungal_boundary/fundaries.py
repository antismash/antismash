# License: GNU Affero General Public License v3 or later
# A copy of GNU AGPL v3 should have been included in this software package in LICENSE.txt.
"""Predict fungal cluster boundaries from motif enrichment in promoters."""


from dataclasses import dataclass, asdict, field
import random
import re
from typing import Any
import numpy as np

from antismash.common.secmet import Gene, Record

from . import utils
from .config import (
    P_VAL_MAX,
    DOWNSTREAM,
    BUFFER_GENES,
    LINEAR_UPSTREAM,
    MAX_COUNTS,
    MAX_GENES,
    MAX_PREDICTED,
    MIN_GENES,
    MIN_GENES_WITH_MOTIFS,
    MIN_LINEAR_SUPPORT,
    MIN_LLR_PER_GENE,
    MODES,
    MOTIF_TYPE_RANGE,
    N_FLANKING,
    N_MOTIFS,
    NUM_TEST_ANCHORS,
    MIN_FLANK,
    SEED,
    TYPES,
    UPSTREAM,
)

# Data Structures
@dataclass
class GeneWithPromoters:  # pylint: disable=too-many-instance-attributes
    """One gene's annotation and promoter sequences for motif searches.

    Both promoters follow the gene's reading direction, with position 0 at the
    upstream end of the extracted sequence. Each offset is the zero-based
    position of the annotated gene start within its corresponding promoter;
    this does not identify the actual transcription start site.
    """
    name: str
    start: int
    end: int
    strand: int
    promoter: str
    tss_offset: int
    linear_promoter: str
    linear_tss_offset: int

    @classmethod
    def from_feature(cls, record: Record, gene: Gene) -> "GeneWithPromoters":
        """Extract promoter sequences and offsets for one annotated gene."""
        start = int(gene.location.start)
        end = int(gene.location.end)
        strand = gene.location.strand
        if strand == 1:
            low = max(0, start - UPSTREAM)
            high = min(len(record.seq), start + DOWNSTREAM)
            promoter = str(record.seq[low:high])
            tss_offset = start - low
        else:
            low = max(0, end - DOWNSTREAM)
            high = min(len(record.seq), end + UPSTREAM)
            promoter = str(record.seq[low:high].reverse_complement())
            tss_offset = high - end
        # Linear motifs are searched closer to the gene than repeats are.
        crop = max(0, tss_offset - LINEAR_UPSTREAM)
        return cls(name=gene.get_name(),
                   start=start, end=end, strand=strand,
                   promoter=promoter, tss_offset=tss_offset,
                   linear_promoter=promoter[crop:],
                   linear_tss_offset=tss_offset - crop)


class Contig:
    """Genes, promoter motif counts, and core anchors for one contig."""

    def __init__(self, record: Record) -> None:
        self.id = record.id
        sorted_genes = sorted(
            record.get_genes(),
            key=lambda feature: (feature.location.start, feature.location.end),
        )
        # composition of the contig, with promoters and motif counts for each gene/promoter
        self.genes: list[GeneWithPromoters] = [
            GeneWithPromoters.from_feature(record, gene_feature)
            for gene_feature in sorted_genes
        ]
        self.gene_motif_counts = [count_motifs(gene) for gene in self.genes]
        self.totals = motif_totals(self.gene_motif_counts)
        self.anchors = core_anchors(record, self.genes)


@dataclass
class GeneRow:  # pylint: disable=too-many-instance-attributes
    """One gene in the evidence table around a chosen motif window."""
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
    mode: str                  # "anchor", "adjacent" or "all"
    type: str                  # "inverted", "direct" or "linear"
    motif: str                 # readable form, e.g. "ACG-N2-CGT"
    motif_id: int
    left: int                  # window start, relative to the anchor
    right: int                 # window end, relative to the anchor
    gene_count: int            # genes in the window
    llr: float                 # enrichment score
    llr_per_gene: float
    sites: int                 # hits in the window (capped counts)
    genes_with_site: int       # genes in the window with the motif
    expected_sites: float      # hits the background predicts
    p_empirical: float         # how often random genes scored as high
    failed_checks: list[str] = field(default_factory=list)


@dataclass
class FundariesReport:
    """The boundary analysis output for one contig."""

    contig: str
    accepted: int
    analysed: int
    results: list[dict[str, Any]]

# Read data and motif counting


def count_motifs(gene: GeneWithPromoters) -> tuple[np.ndarray, np.ndarray]:
    """Return sorted motif IDs and capped counts for this gene's promoters.

    There is a cap on the number of ocurrences in each promoter.
    """

    ids = np.concatenate([utils.repeat_motif_ids(gene.promoter),
                          utils.linear_motif_ids(gene.linear_promoter)])
    counts = np.bincount(ids)
    motifs = np.flatnonzero(counts)
    counts = np.minimum(counts[motifs], MAX_COUNTS)
    return motifs, counts


def motif_totals(gene_motif_counts: list[tuple[np.ndarray, np.ndarray]]) -> np.ndarray:
    """It returns the total of ocurrences of each motif in a contig"""
    totals = np.zeros(N_MOTIFS)
    for motifs, counts in gene_motif_counts:
        totals[motifs] += counts
    return totals


def core_anchors(record: Record, genes: list[GeneWithPromoters]) -> list[tuple[str, int]]:
    """Extract core genes of the record, used as the anchors."""
    positions = {}
    for index, gene in enumerate(genes):
        positions[gene.name] = index
    anchors = []
    for cds in record.get_cds_features():
        if cds.gene_function.name != "CORE":
            continue
        name = cds.get_name()
        index = positions.get(name)
        anchors.append((name, index))
    return anchors


# Scoring functions

def neighbourhood(contig: Contig, index: int, reach: int) -> range:
    """Returns neighbourhood of genes from the given anchor"""
    first = max(0, index - reach)
    last = min(len(contig.genes), index + reach + 1)
    return range(first, last)


def gene_windows(n_local: int, center: int) -> tuple[np.ndarray, np.ndarray]:
    """Returns every window of neighbouring genes that contains the anchor."""
    lefts = []
    rights = []
    for left in range(max(0, center - MAX_GENES + 1), center + 1):
        for right in range(center, min(n_local, left + MAX_GENES)):
            if right - left + 1 >= MIN_GENES:
                lefts.append(left)
                rights.append(right)
    return np.array(lefts), np.array(rights)


def score_windows(contig: Contig, gene_index: int) -> dict:
    """For every motif, find the window of genes around this one where it is
    most enriched compared to the rest of the contig.

    A window scores high when the motif occurs there more often than the
    background rate predicts. The score is a Poisson log-likelihood ratio:
    observed * log(observed / expected) - (observed - expected).

    Returns a dictionary of per-motif arrays (llr, observed, with_site,
    expected, left, right, size) describing each motif's best window,
    with left and right relative to the anchor, plus the anchor's index,
    its row (center) and the neighbourhood count table (counts).
    """
    # Count every motif in each gene near the anchor
    nearby_genes = neighbourhood(contig, gene_index, N_FLANKING)
    center = gene_index - nearby_genes[0]   # anchor's row in the table

    counts = np.zeros((len(nearby_genes), N_MOTIFS))
    for row, gene in enumerate(nearby_genes):
        motifs, hits = contig.gene_motif_counts[gene]
        counts[row, motifs] = hits

    # Calculate the background rate
    genes_outside = len(contig.genes) - len(nearby_genes)
    if genes_outside < 1:
        raise ValueError(f"{contig.id}: leaves no promoters "
                         f"outside the {len(nearby_genes)}-gene neighbourhood "
                         f"to compare against")

    rate = (contig.totals - counts.sum(axis=0)) / genes_outside
    rate = np.maximum(rate, 0.5 / genes_outside)

    # Score possible winners together for each left endpoint.
    lefts, rights = gene_windows(len(nearby_genes), center)
    best = utils.score_motif_windows(counts, rate, lefts, rights, center)
    best_size = best.right - best.left + 1

    return {
        "index": gene_index,
        "center": center,
        "counts": counts,
        "llr": best.llr,
        "observed": best.observed,
        "with_site": best.with_site,
        "expected": best_size * rate,
        "left": best.left - center,
        "right": best.right - center,
        "size": best_size,
    }


def mode_motifs(scored: dict, mode: str) -> np.ndarray:
    """IDs of the motifs a mode may choose from.
    "anchor":   motifs present in the anchor gene's own promoter
    "adjacent": motifs present in the anchor or the genes either side of it
    "all":      every motif
    """
    counts = scored["counts"]
    center = scored["center"]

    if mode == "anchor":
        return np.flatnonzero(counts[center])
    if mode == "adjacent":
        adj_genes = counts[max(0, center - 1):center + 2]
        return np.flatnonzero(adj_genes.sum(axis=0))
    if mode == "all":
        return np.arange(N_MOTIFS)


def highest_llr(scored: dict) -> dict[tuple[str, str], float]:
    """Returns the top score per mode and motif type."""
    all_scores = scored["llr"]
    best = {}
    for mode in MODES:
        mode_ids = mode_motifs(scored, mode)
        for motif_type in TYPES:
            start, stop = MOTIF_TYPE_RANGE[motif_type]
            # Motifs of this mode that belong to this type.
            in_type = (mode_ids >= start) & (mode_ids < stop)
            type_ids = mode_ids[in_type]
            # Highest score among them, or 0 if there are none.
            if len(type_ids) == 0:
                top_score = 0.0
            else:
                top_score = float(all_scores[type_ids].max())
            best[(mode, motif_type)] = top_score
    return best


# calibration against random genes

def calibrate(contig: Contig) -> dict[tuple[str, str], list[float]]:
    """Run the same search at genes to see what chance alone produces.

    Thousands of motif and window combinations are tried per anchor, so a high
    score has to be compared against the high scores reached at arbitrary
    positions on this contig, which carries its own GC content and
    gene density.

    Returns, for each (mode, motif type), the list of best scores reached at
    the random genes, one per gene.
    """
    selected_anchors = range(MIN_FLANK, len(contig.genes) - MIN_FLANK)
    if not selected_anchors:
        raise ValueError(f"{contig.id}: too few genes to calibrate on")
    rng = random.Random(SEED)
    sample_size = min(NUM_TEST_ANCHORS, len(selected_anchors))
    random_genes = rng.sample(selected_anchors, sample_size)

    calibration_runs = {(mode, motif_type): []
                        for mode in MODES for motif_type in TYPES}

    for gene_index in random_genes:
        window_scores = score_windows(contig, gene_index)
        best_per_category = highest_llr(window_scores)

        for category, best_score in best_per_category.items():
            calibration_runs[category].append(best_score)

    return calibration_runs


def empirical_p(llr: float, external_llr: list[float]) -> float:
    """How often by chance reached this score, counting the
    observation itself."""
    times_higher = sum(value >= llr for value in external_llr)
    p_value = (1 + times_higher) / (len(external_llr) + 1)
    return p_value


# Decision module: build, check and rank results
def build_candidates(scored: dict, external_runs: dict) -> list[dict]:
    """Build candidates from all the possible motifs. One candidate per mode
    and enriched motif."""
    candidates = []
    for mode in MODES:
        motifs = mode_motifs(scored, mode)
        for motif in motifs[scored["llr"][motifs] > 0]:
            motif_kind = utils.motif_type(motif)
            llr = float(scored["llr"][motif])
            size = int(scored["size"][motif])
            candidates.append(Candidate(
                mode=mode,
                type=motif_kind,
                motif=utils.motif_int_to_text(motif),
                motif_id=int(motif),
                left=int(scored["left"][motif]),
                right=int(scored["right"][motif]),
                gene_count=size,
                llr=llr,
                llr_per_gene=llr / size,
                sites=int(scored["observed"][motif]),
                genes_with_site=int(scored["with_site"][motif]),
                expected_sites=float(scored["expected"][motif]),
                p_empirical=empirical_p(llr, external_runs[mode, motif_kind]),
            ))
    return candidates


def failed_checks(candidate: Candidate) -> list[str]:
    """The acceptance rules for a candidate to be accepted"""
    failed = []
    if candidate.p_empirical >= P_VAL_MAX:
        failed.append("not_better_than_random_genes")
    if candidate.genes_with_site < MIN_GENES_WITH_MOTIFS:
        failed.append("too_few_promoters_with_motif")
    if candidate.gene_count > MAX_PREDICTED:
        failed.append("window_too_long")
    if candidate.llr_per_gene < MIN_LLR_PER_GENE:
        failed.append("score_per_gene_too_low")
    support = candidate.genes_with_site / candidate.gene_count
    if candidate.type == "linear" and support < MIN_LINEAR_SUPPORT:
        failed.append("linear_motif_missing_from_too_many_promoters")
    return failed


def ranking(candidate: Candidate) -> tuple[bool, float, float, int, str, int]:
    """Sort key for choosing the best candidate: smaller is better"""
    is_linear = candidate.type == "linear"
    p_value = candidate.p_empirical
    negative_score = -candidate.llr
    mode_order = MODES.index(candidate.mode)
    motif_type = candidate.type
    motif_id = candidate.motif_id

    return (is_linear, p_value, negative_score, mode_order,
            motif_type, motif_id)


# Reporting functions

def motif_forms(candidate: Candidate) -> set[str]:
    """Return of a motif in its written from"""
    if candidate.type == "linear":
        word = candidate.motif
        return {word, utils.reverse_complement(word)}
    arm, spacer, other = candidate.motif.split("-")
    gap = "." * int(spacer[1:])
    forms = {arm + gap + other}
    if candidate.type == "direct":
        forms.add(utils.reverse_complement(arm) + gap +
                  utils.reverse_complement(other))
    return forms


def motif_hits(sequence: str, forms: set[str], tss_offset: int) -> list[int]:
    """Returns positions of every occurrence of the motif, relative
    to the gene start."""
    hits = []
    for form in forms:
        matches = re.finditer(f"(?={form})", sequence)
        hits += [match.start() - tss_offset for match in matches]
    return sorted(hits)


def gene_table(contig: Contig, scored: dict,
               candidate: Candidate) -> list[GeneRow]:
    """Returns a summary of where the motif occurs in each gene in its
    window"""
    forms = motif_forms(candidate)
    anchor = scored["index"]
    # The margin stops at the end of the contig.
    first = max(0, anchor + candidate.left - BUFFER_GENES)
    last = min(len(contig.genes) - 1, anchor + candidate.right + BUFFER_GENES)
    rows = []
    for index in range(first, last + 1):
        gene = contig.genes[index]
        if candidate.type == "linear":
            hits = motif_hits(gene.linear_promoter,
                              forms, gene.linear_tss_offset)
        else:
            hits = motif_hits(gene.promoter, forms, gene.tss_offset)

        position = index - anchor
        rows.append(GeneRow(
            name=gene.name,
            position=position,
            in_cluster=candidate.left <= position <= candidate.right,
            start=gene.start,
            end=gene.end,
            strand=gene.strand,
            motif_hits=len(hits),
            hit_positions_from_gene_start=hits,
        ))
    return rows


# Function deciding boundaries of anchor gene

def decide_boundaries(contig: Contig,
                      calibration_runs: dict[tuple[str, str], list[float]],
                      core_gene: str, index: int) -> dict:
    """Decide the cluster boundaries around one core gene.
    Scores every window around the gene,and keeps the enriched
    as candidates and checks against the acceptance
    rules. If any pass, the best one is accepted."""

    scored = score_windows(contig, index)
    candidates = build_candidates(scored, calibration_runs)
    for candidate in candidates:
        candidate.failed_checks = failed_checks(candidate)
    passing = [c for c in candidates if not c.failed_checks]
    result = {"core_gene": core_gene, "contig": contig.id,
              "candidates_enriched": len(candidates),
              "candidates_passing": len(passing)}

    if not passing:
        nearest = min(candidates, key=ranking) if candidates else None
        result["status"] = "rejected"
        result["best_rejected_candidate"] = nearest
        return result

    chosen = min(passing, key=ranking)
    first = index + chosen.left
    last = index + chosen.right
    result["status"] = "accepted"
    result.update(asdict(chosen))
    result["genes"] = [contig.genes[i].name for i in range(first, last + 1)]
    result["gene_table"] = gene_table(contig, scored, chosen)
    return result


def main(record) -> FundariesReport:
    """Main function for running fundaries"""
    # Step 1: read contig, find and count motifs.
    contig = Contig(record)

    # Step 2: Find the best scores random genes reach by chance
    external_analysis = calibrate(contig)

    # Step 3: Score each core gene
    results = []
    for core_gene, index in contig.anchors:
        results.append(decide_boundaries(contig, external_analysis,
                                         core_gene, index))

    # Step 4: Report.
    accepted = []

    for result in results:
        if result["status"] == "accepted":
            accepted.append(result)
    report = FundariesReport(contig=contig.id,
                             accepted=len(accepted),
                             analysed=len(results),
                             results=results)

    return report
