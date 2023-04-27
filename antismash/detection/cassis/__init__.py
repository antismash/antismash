# License: GNU Affero General Public License v3 or later
# A copy of GNU AGPL v3 should have been included in this software package in LICENSE.txt.

"""Implementation of the CASSIS method for the motif-based prediction of SM cluster regions"""

from collections import defaultdict
from functools import reduce
import logging
import operator
import os
import shutil
from typing import Any, Dict, Iterable, List, Optional, Tuple

from Bio.SeqFeature import SeqFeature

from antismash.common import module_results, subprocessing
from antismash.common.secmet import SubRegion, Feature, FeatureLocation, GeneFunction, Gene, Record
from antismash.common.serialiser import feature_to_json, feature_from_json
from antismash.config import ConfigType
from antismash.config.args import ModuleArgs
from antismash.detection import DetectionStage

from .cluster_prediction import get_predictions_for_anchor, ClusterPrediction
from .config import MAX_GAP_LENGTH, MAX_PERCENTAGE
from .pairings import PROMOTER_RANGE
from .promoters import Promoter, CombinedPromoter, get_promoters, DuplicatePromoterError, \
                       write_promoters_to_file

NAME = "cassis"
SHORT_DESCRIPTION = "Detect secondary metabolite gene cluster (motif based)"
DETECTION_STAGE = DetectionStage.AREA_REFINEMENT


class CassisResults(module_results.DetectionResults):
    """ Contains the subregions predicted by cassis """
    def __init__(self, record_id: str) -> None:
        super().__init__(record_id)
        self.subregions: List[SubRegion] = []
        self.promoters: List[Promoter] = []

    def to_json(self) -> Dict[str, Any]:
        subregions = []
        promoters = []

        for cluster in self.subregions:
            subregions.append(feature_to_json(cluster.to_biopython()[0]))

        for promoter in self.promoters:
            promoters.append(promoter.to_json())

        return {"record_id": self.record_id, "subregions": subregions, "promoters": promoters,
                "max_percentage": MAX_PERCENTAGE, "max_gap_length": MAX_GAP_LENGTH}

    @staticmethod
    def from_json(json: Dict[str, Any], record: Record) -> Optional["CassisResults"]:
        # throw away the results if the conditions are different
        if json["record_id"] != record.id:
            logging.debug("Record identifiers don't match, discarding previous results")
            return None
        if json["max_percentage"] != MAX_PERCENTAGE:
            logging.debug("CASSIS commonality threshold changed, discarding previous results")
            return None
        if json["max_gap_length"] != MAX_GAP_LENGTH:
            logging.debug("CASSIS maximum island length changed, discarding previous results")
            return None

        subregions = []
        promoters: List[Promoter] = []
        for cluster in json["subregions"]:
            subregions.append(SubRegion.from_biopython(feature_from_json(cluster)))
        for promoter in json["promoters"]:
            if promoter["type"] == "CombinedPromoter":
                promoters.append(CombinedPromoter.from_json(promoter))
            else:
                promoters.append(Promoter.from_json(promoter))
        results = CassisResults(record.id)
        results.subregions = subregions
        results.promoters = promoters
        return results

    def add_to_record(self, record: Record) -> None:
        store_promoters(self.promoters, record)

    def get_predicted_subregions(self) -> List[SubRegion]:
        # for the purposes of the record, don't return all possible matches,
        # only return non-contained ones
        return filter_subregions(self.subregions)


def get_arguments() -> ModuleArgs:
    """ Build and return arguments. No extra options beyond a switch to enable """
    args = ModuleArgs('Additional analysis', 'cassis')
    args.add_analysis_toggle('--cassis',
                             dest='cassis',
                             action='store_true',
                             default=False,
                             help="Motif based prediction of SM gene cluster regions.")
    return args


def is_enabled(options: ConfigType) -> bool:
    """ Is the module enabled """
    return options.cassis


def check_options(options: ConfigType) -> List[str]:
    """ Make sure the options are sane """
    problems = []
    if options.taxon != "fungi" and is_enabled(options):
        problems.append("CASSIS cluster prediction only works for fungal sequences.")
    return problems


def regenerate_previous_results(previous: Dict[str, Any], record: Record,
                                _options: ConfigType) -> Optional[CassisResults]:
    """ Rebuild the previous run results from a JSON object into this module's
        python results class.

        Arguments:
            previous: the previous results as a dictionary
            record: the Record that was used to generate the previous results
            options: an antismash.Config object
    """
    results = CassisResults.from_json(previous, record)
    if results:
        store_promoters(results.promoters, record)
    return results


def run_on_record(record: Record, results: CassisResults, options: ConfigType) -> CassisResults:
    """ Run this module's analysis section on the given record or use the
        previous results.

        Arguments:
            record: the Record instance to analyse
            results: the previous results as generated by regenerate_previous_results()
            options: an antismash.Config object

        Returns:
            this module's results as a subclass of
                antismash.common.module_results.ModuleResults
    """
    if results:
        logging.debug("Cassis reusing %d cluster(s)", len(results.subregions))
    else:
        results = detect(record, options)
        logging.debug("Cassis detected %d cluster(s)", len(results.subregions))
    results.add_to_record(record)
    return results


def check_prereqs(options: ConfigType) -> List[str]:
    """Check for prerequisites"""
    failure_messages = []

    expected_memesuite_version = (4, 11, 2)
    expected = ".".join(map(str, expected_memesuite_version))
    missing_message = "Failed to locate executable for %r"
    incompatible_message = "Incompatible %s version, expected %s or later but found %s"

    def older_than_expected(version: str) -> bool:
        try:
            if tuple(map(int, version.split("."))) < expected_memesuite_version:
                return True
        except ValueError:
            return True  # since it'll be the same result
        return False

    if "meme" in options.executables:
        meme_version = subprocessing.run_meme_version()
        # later versions of MEME are broken, e.g. 5.5.2 ignores the `-nmotifs` arg
        if meme_version != expected:
            failure_messages.append(f"Incompatible MEME version, expected {expected} but found {meme_version}")
    else:
        failure_messages.append(missing_message % "meme")

    if "fimo" in options.executables:
        fimo_version = subprocessing.run_fimo_version()
        # FIMO changed output format after 4.11.2, but it's simple and handled in the parsers
        if older_than_expected(fimo_version):
            failure_messages.append(incompatible_message % ("FIMO", expected, fimo_version))
    else:
        failure_messages.append(missing_message % "fimo")

    return failure_messages


def detect(record: Record, options: ConfigType) -> CassisResults:
    """Use core genes (anchor genes) from hmmdetect as seeds to detect gene
       cluster regions
    """
    logging.info("Detecting gene cluster regions using CASSIS")

    results = CassisResults(record.id)

    # get core genes from hmmdetect --> necessary CASSIS input, aka "anchor genes"
    anchor_gene_names = get_anchor_gene_names(record)
    logging.debug("Record has %d anchor genes", len(anchor_gene_names))
    if not anchor_gene_names:
        return results

    # filter all genes in record for neighbouring genes with overlapping annotations
    genes = record.get_genes()
    logging.debug("Record has %d features of type 'gene'", len(genes))
    if not genes:
        return results
    candidate_genes, ignored_genes = ignore_overlapping(list(genes))

    # compute promoter sequences/regions --> necessary for motif prediction (MEME and FIMO input)
    try:
        # why these values? see "Wolf et al (2015): CASSIS and SMIPS ..."
        upstream_tss = 1000  # nucleotides upstream TSS
        downstream_tss = 50  # nucleotides downstream TSS
        promoters = get_promoters(record, candidate_genes, upstream_tss, downstream_tss)
        results.promoters = promoters
        write_promoters_to_file(options.output_dir, record.name, promoters)
    except DuplicatePromoterError:
        logging.error("CASSIS discovered an error while working on the promoter sequences, skipping CASSIS analysis")
        return results

    if not promoters:
        logging.debug("CASSIS found zero promoter regions, skipping CASSIS analysis")
        return results
    if len(promoters) < 3:
        logging.debug("Sequence %r yields less than 3 promoter regions, skipping CASSIS analysis", record.name)
        return results

    if len(promoters) < 40:
        logging.debug("Sequence %r yields only %d promoter regions", record.name, len(promoters))
        logging.debug("Cluster detection on small sequences may lead to incomplete cluster predictions")

    predicted_subregions = []
    cluster_predictions = {}  # {anchor gene: cluster predictions}
    for i, anchor in enumerate(anchor_gene_names):
        logging.debug("Detecting cluster around anchor gene %r (%d of %d)", anchor, i + 1, len(anchor_gene_names))
        # get cluster predictions sorted by cluster abundance
        # (most abundant --> "best" prediction)
        predictions = get_predictions_for_anchor(anchor, promoters, record, ignored_genes, options)
        if predictions:
            cluster_predictions[anchor] = predictions
            predicted_subregions.extend(create_subregions(anchor, predictions, record))

    logging.debug("Cleaning up MEME and FIMO output directories")
    cleanup_outdir(anchor_gene_names, cluster_predictions, options)
    results.subregions = predicted_subregions
    return results


def get_anchor_gene_names(record: Record) -> List[str]:
    """ Finds all gene names that have a CDS with secondary metabolite
        annotations.

        Requires that a CDS.get_name() returns the same name of its parent
        Gene.get_name()

        Arguments:
            record: the record to search

        Returns:
            a list of gene names
    """
    anchor_genes = []

    for feature in record.get_cds_features():
        if feature.gene_function == GeneFunction.CORE:
            anchor_genes.append(feature.get_name())

    return anchor_genes


def ignore_overlapping(genes: List[Gene]) -> Tuple[List[Gene], List[Gene]]:
    """Ignore genes with overlapping locations (skip the second gene of an overlapping couple)"""
    ignored = []

    overlap = True
    while overlap:  # check again until no overlap found in the entire (remaining) gene list
        overlap = False
        non_overlapping = [genes[0]]

        for i in range(1, len(genes)):
            if genes[i-1].overlaps_with(genes[i]):
                logging.debug("Ignoring %r (overlapping with %r)",
                              genes[i].get_name(), genes[i-1].get_name())
                ignored.append(genes[i])
                overlap = True
            else:
                non_overlapping.append(genes[i])

        genes = non_overlapping

    if ignored:
        logging.debug("Ignoring %d genes due to overlapping locations", len(ignored))

    return (genes, ignored)


def cleanup_outdir(anchor_gene_names: Iterable[str], cluster_predictions: Dict[str, List[ClusterPrediction]],
                   options: ConfigType) -> None:
    """Delete unnecessary files to free disk space"""
    all_motifs = set()
    for motif in PROMOTER_RANGE:
        all_motifs.add(motif.pairing_string)

    for anchor in anchor_gene_names:
        if anchor in cluster_predictions:
            used_motifs = set()
            for cluster in cluster_predictions[anchor]:
                used_motifs.add(cluster.start.pairing_string)
                used_motifs.add(cluster.end.pairing_string)
            unused_motifs = all_motifs.difference(used_motifs)
            # only remove directories from "unused" motifs (no cluster prediction)
            for directory in unused_motifs:
                shutil.rmtree(os.path.join(options.output_dir, "meme", anchor, directory), ignore_errors=True)
                shutil.rmtree(os.path.join(options.output_dir, "fimo", anchor, directory), ignore_errors=True)
        else:
            # all motifs are "unused" (not a single prediction for this anchor gene)
            # --> remove anchor genes directory, including all motif subdirectories
            shutil.rmtree(os.path.join(options.output_dir, "meme", anchor), ignore_errors=True)
            shutil.rmtree(os.path.join(options.output_dir, "fimo", anchor), ignore_errors=True)


def store_promoters(promoters: Iterable[Promoter], record: Record) -> None:
    """Store information about promoter sequences to a SeqRecord"""
    for promoter in promoters:
        # remember to account for 0-indexed start location
        new_feature = SeqFeature(FeatureLocation(max(0, promoter.start - 1), promoter.end),
                                 type="promoter")
        new_feature.qualifiers = {
            "locus_tag": promoter.get_gene_names(),  # already a list with one or two elements
            "seq": [str(promoter.seq)],
        }

        if isinstance(promoter, CombinedPromoter):
            new_feature.qualifiers["note"] = ["bidirectional promoter"]

        secmet_version = Feature.from_biopython(new_feature)
        secmet_version.created_by_antismash = True

        record.add_feature(secmet_version)


def filter_subregions(subregions: List[SubRegion]) -> List[SubRegion]:
    """ Strips any subregion that is fully contained by another for the same anchor

        Arguments:
            subregions: the subregions to filter

        Returns:
            a sorted list of SubRegions
    """
    if not subregions:
        return subregions

    by_anchor: Dict[str, List[SubRegion]] = defaultdict(list)
    # sort from largest to smallest to avoid complicated replacement logic
    # any sharing an anchor will overlap on that gene anyway
    for sub in sorted(subregions, key=lambda x: x.location.end - x.location.start, reverse=True):
        contained = False
        for other in by_anchor[sub.label]:
            if sub.is_contained_by(other):
                contained = True
                break
        if not contained:
            by_anchor[sub.label].append(sub)
    # flatten the lists and sort back into location order, then anchor
    # mypy doesn't handle this reduce well
    flattened = reduce(operator.concat, by_anchor.values())  # type: ignore
    return sorted(flattened, key=lambda x: (x.location.start, x.location.end, x.label))


def create_subregions(anchor: str, cluster_preds: List[ClusterPrediction],
                      record: Record) -> List[SubRegion]:
    """ Create the predicted subregions """
    subregions: List[SubRegion] = []
    if not cluster_preds:
        return subregions
    for i, cluster in enumerate(cluster_preds):
        # clusters returned by hmmdetect are based on CDS features
        # in contrast, subregions returned by cassis are based on gene features
        # --> hmmdetect derived clusters have exact loctions, like the CDSs have
        # --> cassis derived subregions may have fuzzy locations, like the genes have
        left_name = cluster.start.gene
        right_name = cluster.end.gene
        left = None
        right = None
        for gene in record.get_genes():
            if gene.get_name() == left_name:
                left = gene
            if gene.get_name() == right_name:
                right = gene
            if left and right:
                break
        assert left and right, "boundary genes no longer present in Record"
        new_feature = SeqFeature(
            FeatureLocation(left.location.start, right.location.end), type="subregion")
        new_feature.qualifiers = {
            "aStool": ["cassis"],
            "label": [anchor],
            "abundance": [cluster.start.abundance + cluster.end.abundance],
            "motif_score": [f"{cluster.start.score + cluster.end.score:.1e}"],
            "gene_left": [cluster.start.gene],
            "promoter_left": [cluster.start.promoter],
            "abundance_left": [cluster.start.abundance],
            "motif_left": [cluster.start.pairing_string],
            "motif_score_left": [f"{cluster.start.score:.1e}"],
            "gene_right": [cluster.end.gene],
            "promoter_right": [cluster.end.promoter],
            "abundance_right": [cluster.end.abundance],
            "motif_right": [cluster.end.pairing_string],
            "motif_score_right": [f"{cluster.end.score:.1e}"],
            "genes": [cluster.genes],
            "promoters": [cluster.promoters],
        }

        if i == 0:
            new_feature.qualifiers["note"] = [f"best prediction (most abundant) for anchor gene {anchor}"]
        else:
            new_feature.qualifiers["note"] = [f"alternative prediction ({i}) for anchor gene {anchor}"]

        new_feature = SubRegion.from_biopython(new_feature)
        subregions.append(new_feature)

    return subregions
