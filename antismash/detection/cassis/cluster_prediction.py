# License: GNU Affero General Public License v3 or later
# A copy of GNU AGPL v3 should have been included in this software package in LICENSE.txt.

""" Cluster prediction logic and classes for CASSIS """

import logging
import os
from typing import Any, Dict, List, Optional

from antismash.common.secmet import Record, Gene
from antismash.config import ConfigType

from .config import VERBOSE_DEBUG
from .islands import get_islands, Island
from .motifs import generate_motifs, Motif, filter_meme_results, filter_fimo_results
from .pairings import Pairing
from .promoters import Promoter, get_anchor_promoter_index
from .runners import run_fimo, run_meme

# if alternatives for predictions can contain ignored genes
ALLOW_CLUSTERS_WITH_IGNORED_GENES = False


class ClusterMarker(Pairing):
    """ Marks the start or end position for a cluster by gene and motif.
        Tracks abundance and the motif string (+n/-m) for that position.
    """
    def __init__(self, gene: str, motif: Motif) -> None:
        super().__init__(motif.plus, motif.minus)
        self.gene = str(gene)
        self.abundance = 1
        assert motif.score is not None
        self.score = float(motif.score)
        self.promoter: Optional[str] = None

    def update(self, motif: Motif) -> None:
        """ Uses the given motif's values instead of the old ones, if the given
            motif has a better (lower) score.
        """
        if motif.score is not None and self.score > float(motif.score):
            self.score = float(motif.score)
            self.plus = int(motif.plus)
            self.minus = int(motif.minus)

    def __str__(self) -> str:
        return f"gene {self.gene}; abundance {self.abundance}; motif {self.pairing_string}; score {self.score}"

    def __repr__(self) -> str:
        return f"ClusterMarker({self}, promoter={self.promoter})"

    def __eq__(self, other: Any) -> bool:
        return (isinstance(other, ClusterMarker)
                and all(getattr(self, key) == getattr(other, key) for key in vars(self)))


class ClusterPrediction:
    """ A prediction of a cluster. Contains start and end ClusterMarkers along with
        counts of genes and promoters within that range.
    """
    def __init__(self, start: ClusterMarker, end: ClusterMarker) -> None:
        self.start = start
        self.end = end
        self.genes = 0
        self.promoters = 0

    def __repr__(self) -> str:
        return (
            f"ClusterPrediction(start={self.start}, end={self.end},"
            f" gene_count={self.genes}, promoter_count={self.promoters})"
        )

    def __eq__(self, other: Any) -> bool:
        return (isinstance(other, ClusterPrediction)
                and all(getattr(self, key) == getattr(other, key) for key in vars(self)))


def get_predictions_for_anchor(anchor: str, promoters: List[Promoter], record: Record,
                               ignored_genes: List[Gene], options: ConfigType) -> List[ClusterPrediction]:
    """ Generate the cluster predictions for an anchor gene.

        Arguments:
            anchor: the name of the anchor gene
            promoters: a list of all gene promoters
            record: the record the anchor gene belongs in
            ignore_genes: the genes which are to be ignored
            options: antimash config

        Returns:
            a list of ClusterPrediction instances, one for each predicted
            cluster location
    """
    try:
        anchor_promoter_index = get_anchor_promoter_index(anchor, promoters)
    except ValueError:
        logging.debug("No promoter region for %r, skipping this anchor gene", anchor)
        return []

    # predict motifs with MEME ("de novo")
    meme_dir = os.path.join(options.output_dir, "meme", anchor)
    motifs = generate_motifs(meme_dir, anchor_promoter_index, promoters)
    exit_code = run_meme(meme_dir, options, VERBOSE_DEBUG)
    if exit_code != 0:
        logging.warning("MEME discovered a problem (exit code %d), skipping this anchor gene", exit_code)
        return []

    motifs = filter_meme_results(meme_dir, motifs, anchor)

    if not motifs:
        logging.debug("Could not predict motifs around %r, skipping this anchor gene", anchor)
        return []

    # search motifs with FIMO ("scanning")
    fimo_dir = os.path.join(options.output_dir, "fimo", anchor)
    exit_code = run_fimo(meme_dir, fimo_dir, record, options, VERBOSE_DEBUG)
    if exit_code != 0:
        logging.warning("FIMO discovered a problem (exit code %d), skipping this anchor gene", exit_code)
        return []
    motifs = filter_fimo_results(motifs, fimo_dir, promoters, anchor_promoter_index)

    if not motifs:
        logging.debug("Could not find motif occurrences for %r, skipping this anchor gene", anchor)
        return []

    # TODO: as4, SiTaR (http://bioinformatics.oxfordjournals.org/content/27/20/2806):
    # Alternative to MEME and FIMO. Part of the original CASSIS implementation.
    # No motif prediction (no MEME). Motif search with SiTaR (instead if FIMO).
    # Have to provide a file in FASTA format with binding site sequences of at least one transcription factor.
    # Will result in binding sites per promoter (like FIMO) --> find islands
    #
    # implement: YES? NO?

    # find islands of binding sites around anchor gene
    islands = get_islands(anchor_promoter_index, motifs, promoters)
    logging.debug("%d possible cluster predictions for %r", len(islands), anchor)
    predictions = create_predictions(islands)
    return check_cluster_predictions(predictions, record, promoters, ignored_genes)


def check_cluster_predictions(cluster_predictions: List[ClusterPrediction],
                              record: Record, promoters: List[Promoter],
                              ignored_genes: List[Gene]) -> List[ClusterPrediction]:
    """Get some more info about each cluster prediction and check if it is sane"""
    checked_predictions = []
    for cp, prediction in enumerate(cluster_predictions):
        sane = True

        # find indices of first and last GENE of the cluster prediction in all genes
        all_gene_names = [gene.get_name() for gene in sorted(record.get_genes(), key=lambda x: x.location.start)]

        if not all_gene_names:
            continue

        start_index_genes = None
        end_index_genes = None
        for i, gene_name in enumerate(all_gene_names):
            if not start_index_genes and prediction.start.gene == gene_name:
                start_index_genes = i
            if not end_index_genes and prediction.end.gene == gene_name:
                end_index_genes = i
            if start_index_genes and end_index_genes:
                break
        assert start_index_genes is not None
        assert end_index_genes is not None

        # find indices of first and last PROMOTER of the cluster prediction in all promoters
        start_index_promoters = None
        end_index_promoters = None
        for i, promoter in enumerate(promoters):
            if not start_index_promoters and prediction.start.gene in promoter.get_gene_names():
                start_index_promoters = i
            if not end_index_promoters and prediction.end.gene in promoter.get_gene_names():
                end_index_promoters = i
            if start_index_promoters and end_index_promoters:
                break
        assert start_index_promoters is not None
        assert end_index_promoters is not None

        prediction.start.promoter = promoters[start_index_promoters].get_id()
        prediction.end.promoter = promoters[end_index_promoters].get_id()
        prediction.genes = end_index_genes - start_index_genes + 1
        prediction.promoters = end_index_promoters - start_index_promoters + 1
        if VERBOSE_DEBUG:
            if cp == 0:
                logging.debug("Best prediction (most abundant): %r -- %r",
                              prediction.start.gene, prediction.end.gene)
            else:
                logging.debug("Alternative prediction (%s): %r -- %r",
                              cp, prediction.start.gene, prediction.end.gene)

        # warn if cluster prediction right at or next to record (~ contig) border
        if start_index_genes < 10:
            if VERBOSE_DEBUG:
                logging.debug("Upstream cluster border located at or next to sequence record border,"
                              " prediction could have been truncated by record border")
            sane = False
        if end_index_genes > len(all_gene_names) - 10:
            if VERBOSE_DEBUG:
                logging.debug("Downstream cluster border located at or next to sequence record border,"
                              " prediction could have been truncated by record border")
            sane = False

        # warn if cluster prediction too short (includes less than 3 genes)
        if prediction.genes < 3:
            if VERBOSE_DEBUG:
                logging.debug("Cluster is very short (less than 3 genes). Prediction may be questionable.")
            sane = False

        # warn if ignored gene (overlapping with anthor gene, see ignore_overlapping())
        # would have been part of the cluster
        for ignored_gene in ignored_genes:
            if ignored_gene.get_name() in all_gene_names[start_index_genes: end_index_genes + 1]:
                if VERBOSE_DEBUG:
                    logging.debug("Gene %r is part of the predicted cluster,"
                                  " but it is overlapping with another gene and was ignored", ignored_gene)
                    logging.debug("Gene %r could have affected the cluster prediction", ignored_gene)
                if not ALLOW_CLUSTERS_WITH_IGNORED_GENES:
                    sane = False
                break

        checked_predictions.append(prediction)

        if sane:  # TODO: this seems wrong
            break
    return checked_predictions


def create_predictions(islands: List[Island]) -> List[ClusterPrediction]:
    """ Sort upstream (start) and downstream (end) borders of islands by abundance
        and create ClusterPredictions from the results.
    """
    # count border abundance
    # start/end are treated independently!
    starts: Dict[str, ClusterMarker] = {}
    ends: Dict[str, ClusterMarker] = {}
    for island in islands:
        if island.start.gene_name in starts:
            scorer = starts[island.start.gene_name]
            scorer.abundance += 1
            scorer.update(island.motif)
        else:
            starts[island.start.gene_name] = ClusterMarker(island.start.gene_name, island.motif)

        if island.end.get_gene_names()[-1] in ends:
            scorer = ends[island.end.get_gene_names()[-1]]
            scorer.abundance += 1
            scorer.update(island.motif)
        else:
            ends[island.end.get_gene_names()[-1]] = ClusterMarker(island.end.get_gene_names()[-1], island.motif)

    # compute sum of start and end abundance, remove duplicates, sort descending
    abundances_sum_sorted = sorted(set(s.abundance + e.abundance
                                       for s in starts.values() for e in ends.values()), reverse=True)
    # compute sum of start and end motif score, remove duplicates, sort ascending
    scores_sum_sorted = sorted(set(s.score + e.score
                                   for s in starts.values() for e in ends.values()))
    # sort by value (=abundance) of start, descending
    starts_sorted = sorted(starts, key=lambda x: starts[x].abundance, reverse=True)
    # sort by value (=abundance) of end, descending
    ends_sorted = sorted(ends, key=lambda x: ends[x].abundance, reverse=True)

    clusters = []
    for abundance in abundances_sum_sorted:
        # list from highest (best) to lowest (worst) abundance
        for score in scores_sum_sorted:
            # list from lowest (best) to highest (worst) motif score/e-value
            for start in starts_sorted:
                for end in ends_sorted:
                    if (starts[start].abundance + ends[end].abundance == abundance
                            and starts[start].score + ends[end].score == score):
                        start_marker = starts[start]
                        end_marker = ends[end]
                        clusters.append(ClusterPrediction(start_marker, end_marker))
                        if VERBOSE_DEBUG:
                            logging.debug("Upstream border:   %s", start_marker)
                            logging.debug("Downstream border: %s", end_marker)
                            logging.debug("Total abundance %s, total score %.1e", abundance, score)
    return clusters
