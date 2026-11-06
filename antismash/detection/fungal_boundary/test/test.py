# License: GNU Affero General Public License v3 or later
# A copy of GNU AGPL v3 should have been included in this software package in
# LICENSE.txt.

"""Tests for fundaries_v7.

    Run from the repository root with:

        python -m pytest antismash/detection/fungal_boundary/test
"""

import math
import unittest

import numpy as np

from antismash.common.secmet import GeneFunction
from antismash.common.secmet.test.helpers import (
    DummyCDS,
    DummyRecord,
)

from .. import config
from .. import fundaries
from .. import utils

GENE_LENGTH_IN_BASES = 600
GENE_START_SPACING_IN_BASES = 2000


def gene_start(index: int) -> int:
    """ Returns the start coordinate of the gene at the given index

        Arguments:
            index: the gene's position on the contig

        Returns:
            the coordinate at which that gene begins
    """
    return config.UPSTREAM + index * GENE_START_SPACING_IN_BASES


def build_dummy_contig(num_genes: int, motif: str = "",
                       motif_gene_indices: list[int] | tuple[int, ...] | range = (), *,
                       cores: list[int] | tuple[int, ...] | range = (),
                       upstream_dist_to_motif: int = 800, seed: int = 0) -> DummyRecord:
    """ Builds a record of evenly spaced genes on random DNA

        Arguments:
            num_genes: number of genes in dummy record contig.
            motif: DNA string of motif added.
            motif_gene_indices: the indices of the genes to plant the motif in.
            cores: the indices of the genes to mark as biosynthetic.
            upstream_dist_to_motif: the distance from gene start to the first motif
            seed: seed for reproducablity

        Returns:
            DummyRecord
    """
    rng = np.random.default_rng(seed)
    sequence = rng.choice(list("ACGT"), size=gene_start(num_genes) + GENE_LENGTH_IN_BASES).tolist()

    for index in motif_gene_indices:
        promoter_start = gene_start(index) - upstream_dist_to_motif
        for offset in (0, 20, 40):
            motif_start = promoter_start + offset
            # inserts the motif
            sequence[motif_start:motif_start + len(motif)] = list(motif)

    record = DummyRecord(seq="".join(sequence))
    for index in range(num_genes):
        start = gene_start(index)
        cds = DummyCDS(start=start, end=start + GENE_LENGTH_IN_BASES,
                       strand=1, locus_tag=f"g{index}")
        if index in cores:
            cds.gene_functions.add(GeneFunction.CORE, "test", "planted anchor", product="test")
        record.add_cds_feature(cds)
    return record


def make_genes(counts_by_gene: list[dict[int, int]], anchor_indices: tuple[int, ...] = ()
               ) -> dict[str, fundaries.GeneWithMotifs]:
    """Build genes from explicit motif counts for scoring tests.

        Arguments:
            counts_by_gene: each gene's motif IDs mapped to their full counts
            anchor_indices: the indices of genes to mark as biosynthetic anchors

        Returns:
            genes in contig order, keyed by gene name
    """
    genes = {}
    for index, motif_counts in enumerate(counts_by_gene):
        motif_ids = sorted(motif_counts)
        counts = np.array([motif_counts[motif_id] for motif_id in motif_ids], dtype=np.int64)
        start = gene_start(index)
        gene = fundaries.GeneWithMotifs(
            name=f"gene_{index}",
            start=start,
            end=start + GENE_LENGTH_IN_BASES,
            strand=1,
            promoter="",
            promoter_to_cds_start_dist=0,
            linear_promoter="",
            linear_promoter_to_cds_start_dist=0,
            motif_ids=np.array(motif_ids, dtype=np.int64),
            counts=np.minimum(counts, config.MAX_COUNTS),
            uncapped_counts=counts,
            is_anchor=index in anchor_indices,
        )
        genes[gene.name] = gene
    return genes


def brute_force_motif_llr(window: list[fundaries.GeneWithMotifs],
                          background_frequency: float, motif: int, centre: int) -> float:
    """ Returns the best eligible LLR for one motif around the anchor

        Arguments:
            window: the genes surrounding the anchor
            background_frequency: the motif's background frequency within its type
            motif: the motif ID to score
            centre: the index of the anchor within the window

        Returns:
            the highest LLR of an eligible window for this motif,
            or 0.0 when no window is eligible
    """
    kind = utils.get_motif_type(motif)
    # per gene: hits of this motif, and hits of every motif of its type
    hits = [int(gene.counts[gene.motif_ids == motif].sum()) for gene in window]
    totals = [
        sum(int(count) for identifier, count in zip(gene.motif_ids, gene.counts)
            if utils.get_motif_type(int(identifier)) == kind)
        for gene in window
    ]

    best = 0.0
    # every window that contains the anchor, within the size limits
    for first in range(max(0, centre - config.MAX_GENES_IN_BGC + 1), centre + 1):
        last = min(len(window), first + config.MAX_GENES_IN_BGC)
        for stop in range(max(centre + 1, first + config.MIN_GENES_WITH_MOTIF), last + 1):
            window_hits = hits[first:stop]
            observed = sum(window_hits)
            expected = sum(totals[first:stop]) * background_frequency
            genes_with_hits = sum(1 for hit in window_hits if hit > 0)
            if genes_with_hits < config.MIN_GENES_WITH_MOTIF or observed <= expected:
                continue
            best = max(best, observed * np.log(observed / expected) - (observed - expected))
    return best


def brute_force_llr(genes: dict[str, fundaries.GeneWithMotifs],
                    background: dict[int, float], anchor: int) -> float:
    """ Returns the best eligible LLR across all motif types for one anchor

        Every window total is summed from scratch, so this agrees with
        find_enriched_motifs only if that function's
        cumulative sums are right.

        Arguments:
            genes: the contig's genes
            background: each motif id's background frequency
            anchor: the index of the anchor gene to score windows around

        Returns:
            the highest LLR of any eligible window containing that anchor,
            or 0.0 when no window is eligible
    """
    gene_list = list(genes.values())
    window = gene_list[max(0, anchor - config.MAX_GENES_IN_BGC):
                       anchor + config.MAX_GENES_IN_BGC + 1]
    centre = anchor - max(0, anchor - config.MAX_GENES_IN_BGC)
    motif_ids = {int(identifier) for gene in window for identifier in gene.motif_ids}
    return max(
        (brute_force_motif_llr(window, background.get(motif, 0.0), motif, centre) for motif in motif_ids),
        default=0.0,
    )


class TestFundaries(unittest.TestCase):
    """ Tests for the fungal boundary pipeline: boundary detection for each
        motif type, the LLR calculation, the length of a reported window,
        the acting anchor and the reproducibility.
    """

    def assert_cluster(self, record: DummyRecord, motif_gene_indices: range,
                       motif_type: str, pattern: str) -> None:
        """ Asserts that a record yields exactly one expected prediction

            Arguments:
                record: the record to run the pipeline over
                motif_gene_indices: the indices of the genes the motif was planted in
                motif_type: the motif class the prediction should name
                pattern: the motif's DNA string, as the report formats it

            Returns:
                None
        """
        predictions = fundaries.main(record).predictions
        self.assertEqual(len(predictions), 1)
        self.assertEqual(predictions[0]["motif_type"], motif_type)
        self.assertEqual(predictions[0]["motif_string_dna"], pattern)
        self.assertEqual(predictions[0]["start"], gene_start(motif_gene_indices[0]))
        self.assertEqual(predictions[0]["end"],
                         gene_start(motif_gene_indices[-1]) + GENE_LENGTH_IN_BASES)

    def test_inverted_repeat_cluster_is_detected(self) -> None:
        """ An inverted repeat cluster is found with the right boundaries

            ACG-N5-CGT: inverted repeat motif to find.
        """
        motif_gene_indices = range(600, 612)
        record = build_dummy_contig(1200, "ACGTTTTTCGT", motif_gene_indices, cores=[606])
        self.assert_cluster(record, motif_gene_indices, "inverted", "ACG-N5-CGT")

    def test_direct_repeat_cluster_is_detected(self) -> None:
        """ A direct repeat cluster is found with the right boundaries

            ACG-N4-ACG: direct repeat motif to find.
        """
        motif_gene_indices = range(600, 612)
        record = build_dummy_contig(1200, "ACGAAAAACG", motif_gene_indices, cores=[606])
        self.assert_cluster(record, motif_gene_indices, "direct", "ACG-N4-ACG")

    def test_linear_motif_cluster_is_detected(self) -> None:
        """ A linear motif cluster is found with the right boundaries

            ACGGTA: linear motif to find.
        """
        motif_gene_indices = range(600, 612)
        record = build_dummy_contig(1200, "ACGGTA", motif_gene_indices, cores=[606],
                                    upstream_dist_to_motif=250)
        self.assert_cluster(record, motif_gene_indices, "linear", "ACGGTA")

    def test_llr_calculation(self) -> None:
        """Check direct sums and cumulative sums for each motif type."""
        expected_llr = 9 * math.log(9 / 0.12) - (9 - 0.12)
        for motif in (1, config.DIRECT_MOTIF_ID_START, config.LINEAR_MOTIF_ID_START):
            with self.subTest(motif=motif):
                # The middle three genes have nine hits out of twelve counts
                # of this type. The flanks lower the score of longer windows.
                flanks = {motif + 1: 3}
                counts = {motif: 3, motif + 1: 1}
                genes = make_genes([flanks, counts, counts, counts, flanks],
                                   anchor_indices=(2,))
                background = {motif: 0.01, motif + 1: 0.99}

                result = fundaries.find_enriched_motifs(
                    genes, background)["gene_2"]

                self.assertAlmostEqual(brute_force_llr(genes, background, 2), expected_llr)
                self.assertIsNotNone(result)
                assert result is not None
                self.assertAlmostEqual(result.llr, expected_llr)
                self.assertEqual(result.motif_id, motif)
                self.assertEqual(result.first_gene_offset, -1)
                self.assertEqual(result.gene_count, 3)

    def test_llr_with_competing_motif_types(self) -> None:
        """Rank repeat and linear motifs by score, preferring repeats on ties."""
        linear_motif = config.LINEAR_MOTIF_ID_START
        expected_llr = 9 * math.log(9 / 0.12) - (9 - 0.12)
        for repeat_motif in (1, config.DIRECT_MOTIF_ID_START):
            examples = (
                (1, 3, 0.2, 0.01, linear_motif),
                (3, 1, 0.01, 0.2, repeat_motif),
                (3, 3, 0.01, 0.01, repeat_motif),
            )
            for repeat_hits, linear_hits, repeat_freq, linear_freq, expected_motif in examples:
                with self.subTest(repeat_motif=repeat_motif, repeat_freq=repeat_freq,
                                  linear_freq=linear_freq):
                    counts = {
                        repeat_motif: repeat_hits, repeat_motif + 1: 4 - repeat_hits,
                        linear_motif: linear_hits, linear_motif + 1: 4 - linear_hits,
                    }
                    flanks = {repeat_motif + 1: 3, linear_motif + 1: 3}
                    genes = make_genes([flanks, counts, counts, counts, flanks],
                                       anchor_indices=(2,))
                    background = {
                        repeat_motif: repeat_freq, repeat_motif + 1: 1 - repeat_freq,
                        linear_motif: linear_freq, linear_motif + 1: 1 - linear_freq,
                    }

                    result = fundaries.find_enriched_motifs(
                        genes, background)["gene_2"]

                    self.assertAlmostEqual(brute_force_llr(genes, background, 2), expected_llr)
                    self.assertIsNotNone(result)
                    assert result is not None
                    self.assertAlmostEqual(result.llr, expected_llr)
                    self.assertEqual(result.motif_id, expected_motif)
                    self.assertEqual(result.first_gene_offset, -1)
                    self.assertEqual(result.gene_count, 3)

    def test_window_length_matches_planted_length(self) -> None:
        """ Tests if the reported window of a run is the expected one.
            In this case, we expect the exac genes with the motif
            planted.

            LLR penalises motif free genes, so the abrupt stop of
            motifs should give the exact result we expect.
        """
        for motif_gene_count in (8, 16):
            record = build_dummy_contig(1200, "ACGTTTTTCGT",
                                        range(600, 600 + motif_gene_count),
                                        cores=[600 + motif_gene_count // 2])
            predictions = fundaries.main(record).predictions
            with self.subTest(motif_gene_count=motif_gene_count):
                self.assertEqual(len(predictions), 1)
                self.assertEqual(len(predictions[0]["gene_table"]), motif_gene_count)

    def test_contig_has_enough_acting_anchors(self) -> None:
        """ Tests if a contig supplies enough acting anchors to
            reach config.MIN_RANDOM_ANCHORS
        """
        cores = range(200, 400, 40)
        record = build_dummy_contig(1200, cores=cores)
        genes = {}
        for cds in record.get_cds_features():
            gene = fundaries.GeneWithMotifs.from_feature(record, cds)
            genes[gene.name] = gene
        results: dict[str, fundaries.EnrichedMotif | None] = {f"g{anchor}": fundaries.EnrichedMotif(
                       anchor_id=f"g{anchor}", motif_id=0, llr=1.0, mode=fundaries.SelectionMode.ALL,
                       first_gene_offset=-15, gene_count=30)
                   for anchor in cores}

        fundaries.select_random_acting_anchors(genes, results)

        selected = sum(gene.used_as_random_acting_anchor
                       for gene in genes.values())
        #needed = next(n for n in range(1, 100)
        #              if 1 / (n + 1) < config.P_VAL_MAX)
        self.assertGreaterEqual(
            selected, config.MIN_RANDOM_ANCHORS,
            f"only {selected} acting anchors; {config.MIN_RANDOM_ANCHORS}"
            "are needed")

    def test_random_anchors_are_reproducible(self) -> None:
        """ Test if with the same seed, the same results are
            found.
        """
        chosen = []
        for _ in range(2):
            # 600 genes, so the contig is big enough that some are picked
            record = build_dummy_contig(600, cores=[300])
            genes = {}
            for cds in record.get_cds_features():
                gene = fundaries.GeneWithMotifs.from_feature(record, cds)
                genes[gene.name] = gene
            fundaries.select_random_acting_anchors(genes, {})
            chosen.append([name for name, gene in genes.items()
                           if gene.used_as_random_acting_anchor])

        self.assertEqual(chosen[0], chosen[1])
