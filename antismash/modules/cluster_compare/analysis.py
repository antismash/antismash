# License: GNU Affero General Public License v3 or later
# A copy of GNU AGPL v3 should have been included in this software package in LICENSE.txt.

""" The core analysis functions for regions as a whole.

    Functions for complicated metrics are split into separate files.
"""

from collections import defaultdict
import logging
import math
import os
from tempfile import NamedTemporaryFile
from typing import (
    Any,
    Dict,
    Iterable,
    List,
    Set,
    Tuple,
    TypeVar,
)

from antismash.common import fasta, subprocessing
from antismash.common.secmet import (
    CDSFeature,
    Record,
)
from antismash.common.secmet.features import CDSCollection, Region, Protocluster

from .components import calculate_component_score, gather_query_components, Components
from .data_structures import (
    load_data,
    DBConfig,
    Hit,
    HitsByReference,
    Mode,
    ReferenceProtocluster,
    ReferenceRecord,
    ReferenceRegion,
    ReferenceScorer,
    ScoresByRegion,
)
from .ordering import calculate_order_score
from .results import (
    DatabaseResults,
    RegionToRegionScores,
    ProtoToRegionScores,
    ProtoToProtoScores,
    VariantResults,
)

HitsByCDS = Dict[CDSFeature, Dict[str, List[Hit]]]
HitsByReferenceName = Dict[str, Dict[str, List[Hit]]]
Ranking = List[ReferenceScorer]
T = TypeVar("T", bound="Dict[ReferenceRegion, Any]")

RANKING_EXPONENT = 3
NORMALISER = math.e**RANKING_EXPONENT - 1
METRIC_CUTOFF = 0.05
MAX_RESULTS = 50


def filter_by_query_area(area: CDSCollection, hits_by_reference: HitsByReference) -> HitsByReference:
    """ Filters a result set for a reference area to those within a specific
        CDSCollection.

        Arguments:
            area: the CDSCollection for filtering
            hits_by_reference: a dictionary mapping ReferenceRecord to
                                a dictionary mapping reference CDS name to Hit

        Returns:
            a dictionary mapping ReferenceRecord to
                a dictionary mapping reference CDS name to Hit
    """
    hits_for_area: HitsByReference = defaultdict(lambda: defaultdict(list))
    for ref_area, hits_by_ref_area in hits_by_reference.items():
        for ref_id, hits in hits_by_ref_area.items():
            for hit in hits:
                if hit.cds in area:
                    hits_for_area[ref_area][ref_id].append(hit)
    return hits_for_area


def filter_by_reference_protocluster(area: ReferenceProtocluster, hits_by_reference: HitsByReference
                                     ) -> HitsByReference:
    """ Filters a result set for a ReferenceRegion to those within a specific
        ReferenceProtocluster.

        Arguments:
            area: the ReferenceProtocluster for filtering
            hits_by_reference: a dictionary mapping ReferenceRecord to
                                a dictionary mapping reference CDS name to Hit

        Returns:
            a dictionary mapping ReferenceRecord to
                a dictionary mapping reference CDS name to Hit
    """
    result: HitsByReference = defaultdict(dict)
    for ref_area, hits_by_ref_area in hits_by_reference.items():
        for ref_cds, hits in hits_by_ref_area.items():
            if ref_cds in area.cdses:
                result[ref_area][ref_cds] = hits
    return result


def score_query_area(query_area: CDSCollection, hits_by_reference: HitsByReference,
                     query_components: Components, mode: Mode) -> Ranking:
    """ Scores a query area against reference areas

        Arguments:
            query_area: the CDSCollection instance to use
            hits_by_reference: a dictionary mapping ReferenceRecord to
                                a dictionary mapping reference CDS name to Hit
            query_components: a dictionary mapping the region and each contained protocluster to
                                a Components instance with the relevant data
            mode: the Mode in which to run the analysis

        Returns:
            a list of ReferenceScorer instances, ranked by decreasing score
    """
    local_hits = filter_by_query_area(query_area, hits_by_reference)

    if not local_hits:
        return []

    best_hits = {ref: trim_to_best_hit(hits) for ref, hits in local_hits.items()}

    query_cdses = query_area.cds_children
    scores: Ranking = []
    for reference_area, hits in best_hits.items():
        identity = calculate_identity_score(hits.values(), len(hits))
        if identity < METRIC_CUTOFF:
            continue
        order = 1.  # can never be less if there's only one thing to hit
        if len(reference_area.cdses) > 1:
            order = calculate_order_score(query_cdses, hits, reference_area.cdses)
            if order < METRIC_CUTOFF:
                continue
        component = calculate_component_score(query_components, reference_area, mode)

        if component is not None and component < METRIC_CUTOFF:
            continue
        scores.append(ReferenceScorer(hits, reference_area, identity, order, component, mode))
    return sorted(scores, reverse=True)


def calculate_protocluster_ranking(scorer: ReferenceScorer) -> float:
    """ Calculates the score for a reference protocluster for use in aggregate
        ranking of reference regions.

        Arguments:
            scorer: the ReferenceScorer result to calculate with

        Returns:
            a float, bounded between 0 and 1, inclusive
    """
    return (math.e**(RANKING_EXPONENT*scorer.final_score) - 1) / NORMALISER


def score_as_protoclusters(label: str, region: Region, hits_by_reference: HitsByReference,
                           query_components: Dict[CDSCollection, Components], mode: Mode
                           ) -> VariantResults:
    """ Performs a protocluster vs reference region comparison

        Arguments:
            label: the name to attach to the results
            region: the query Region
            hits_by_reference: a dictionary mapping ReferenceRecord to
                                a dictionary mapping reference CDS name to Hit
            query_components: a dictionary mapping the region and each contained protocluster to
                                a Components instance with the relevant data
            mode: the Mode in which to run the analysis

        Returns:
            a VariantResults instance
    """
    local_hits = filter_by_query_area(region, hits_by_reference)

    total_scores: Dict[ReferenceRegion, float] = defaultdict(float)

    scores: Dict[int, Dict[ReferenceRegion, ReferenceScorer]] = defaultdict(dict)
    for protocluster in region.get_unique_protoclusters():
        for scorer in score_query_area(protocluster, local_hits, query_components[protocluster], mode):
            total_scores[scorer.reference] += calculate_protocluster_ranking(scorer)
            scores[protocluster.get_protocluster_number()][scorer.reference] = scorer

    ranking = sorted(total_scores.items(), key=lambda x: x[1], reverse=True)
    ranking, scores, best_hits = apply_limits_to_rankings(ranking, scores, local_hits)
    return VariantResults(label, ranking, ProtoToRegionScores(scores), best_hits)


def score_as_region(label: str, region: Region, hits_by_reference: HitsByReference,
                    query_components: Dict[CDSCollection, Components], mode: Mode
                    ) -> VariantResults:
    """ Performs a region vs region comparison

        Arguments:
            label: the name to attach to the results
            region: the query Region
            hits_by_reference: a dictionary mapping ReferenceRecord to
                                a dictionary mapping reference CDS name to Hit
            query_components: a dictionary mapping the region and each contained protocluster to
                                a Components instance with the relevant data
            mode: the Mode in which to run the analysis

        Returns:
            a VariantResults instance
    """
    local_hits = filter_by_query_area(region, hits_by_reference)
    ranking = score_query_area(region, local_hits, query_components[region], mode)[:MAX_RESULTS]
    region_ranking = sorted(((scorer.reference, scorer.final_score) for scorer in ranking),
                            key=lambda x: x[1], reverse=True)
    regions = {region for region, _ in region_ranking}
    best_hits = {ref: trim_to_best_hit(hits) for ref, hits in local_hits.items() if ref in regions}
    return VariantResults(label + str(mode), region_ranking, RegionToRegionScores(ranking), best_hits)


def score_against_protoclusters(label: str, region: Region, hits_by_reference: HitsByReference,
                                query_components: Dict[CDSCollection, Components], mode: Mode
                                ) -> VariantResults:
    """ Performs a protocluster vs protocluster comparison

        Arguments:
            label: the name to attach to the results
            region: the query Region
            hits_by_reference: a dictionary mapping ReferenceRecord to
                                a dictionary mapping reference CDS name to Hit
            query_components: a dictionary mapping the region and each contained protocluster to
                                a Components instance with the relevant data
            mode: the Mode in which to run the analysis

        Returns:
            a VariantResults instance
    """
    score_matrix: Dict[int, Dict[ReferenceRegion,
                                 Dict[ReferenceProtocluster, ReferenceScorer]
                                 ]] = defaultdict(lambda: defaultdict(dict))
    reference_best_scores: Dict[Protocluster, Dict[ReferenceRegion, float]] = defaultdict(lambda: defaultdict(float))
    local_hits = filter_by_query_area(region, hits_by_reference)
    for ref_region in local_hits:
        hits_for_ref_region = {ref_region: local_hits[ref_region]}
        for ref_protocluster in ref_region.protoclusters:
            hits = filter_by_reference_protocluster(ref_protocluster, hits_for_ref_region)
            for protocluster in region.get_unique_protoclusters():
                for scorer in score_query_area(protocluster, hits, query_components[protocluster], mode):
                    score = max(scorer.final_score, reference_best_scores[protocluster][ref_region])
                    reference_best_scores[protocluster][ref_region] = score
                    score_matrix[protocluster.get_protocluster_number()][ref_region][ref_protocluster] = scorer

    reference_total_scores: Dict[ReferenceRegion, float] = defaultdict(float)
    for ref_region_to_score in reference_best_scores.values():
        for ref_region, score in ref_region_to_score.items():
            reference_total_scores[ref_region] += score

    region_ranking = sorted(reference_total_scores.items(), key=lambda x: x[1], reverse=True)

    region_ranking, score_matrix, best_hits = apply_limits_to_rankings(region_ranking, score_matrix, local_hits)
    return VariantResults(label, region_ranking, ProtoToProtoScores(score_matrix), best_hits)


def apply_limits_to_rankings(region_ranking: ScoresByRegion, score_matrix: Dict[int, T],
                             hits: HitsByReference, max_reference_regions: int = MAX_RESULTS,
                             ) -> Tuple[ScoresByRegion, Dict[int, T],
                                        Dict[ReferenceRegion, Dict[str, Hit]]]:
    """ Returns a filtered version of the input matrix, where only the top
        hits are kept.

        Arguments:
            region_ranking: a list of tuples of ReferenceRegion and float score
            score_matrix: a dictionary mapping query protocluster number to
                            a dictionary mapping ReferenceRegion to
                                some subset of results
            hits: the hits for the relevant Region and ReferenceRegion
            max_reference_regions: the maximum number of reference regions to
                                   include in the ranking
        Returns:
            a tuple of
                the first values of region_ranking, and
                the filtered copy of score_matrix
                a dictionary mapping each ReferenceRegion to
                    a dictionary of each CDS to it's singular best hit
    """
    region_ranking = region_ranking[:max_reference_regions]
    regions = {rank[0] for rank in region_ranking}
    limited_matrix: Dict[int, Any] = {}
    for key, pairs_by_ref_region in score_matrix.items():
        limited_matrix[key] = {ref: value for ref, value in pairs_by_ref_region.items() if ref in regions}
    best_hits = {ref: trim_to_best_hit(hits) for ref, hits in hits.items() if ref in regions}
    return region_ranking, limited_matrix, best_hits


def run(record: Record, config: DBConfig) -> DatabaseResults:
    """ Runs the analysis on the given record with the given reference dataset.

        Arguments:
            record: the Record instance to analyse
            data_dir: the directory containing the data to analyse with

        Returns:
            a DatabaseResults instance
    """
    if not record.get_regions():
        return DatabaseResults(config.name, config.url, config.description, {})
    logging.debug("Loading comparison database: %s", config.name)
    references = load_data(os.path.join(config.path, "data.json"))
    logging.debug("Comparison database loaded")

    # find all the hits between this record and reference databases
    _, hits_by_name = find_diamond_matches(record, os.path.join(config.path, "proteins.dmnd"))
    hits = convert_to_references(hits_by_name, references)

    # split into regions and/or protoclusters as required
    results_per_region: Dict[int, Dict[str, VariantResults]] = {}
    for region in record.get_regions():
        # gather query components once and reuse throughout
        query_components: Dict[CDSCollection, Components] = {}
        query_components[region] = gather_query_components(region.cds_children)
        for protocluster in region.get_unique_protoclusters():
            query_components[protocluster] = gather_query_components(protocluster.cds_children)

        # run the requested comparisons
        region_results = {}
        for comparison in config.comparisons:
            label = f"{comparison}_{config.mode}"
            if comparison == "ProtoToProto":
                result = score_against_protoclusters(label, region, hits, query_components, config.mode)
            elif comparison == "ProtoToRegion":
                result = score_as_protoclusters(label, region, hits, query_components, config.mode)
            else:
                result = score_as_region(label, region, hits, query_components, config.mode)
            region_results[label] = result

        results_per_region[region.get_region_number()] = region_results

    return DatabaseResults(config.name, config.url, config.description, results_per_region)


def calculate_identity_score(hits: Iterable[Hit], count: int) -> float:
    """ Calculates a score between 0 and 1, inclusive, of the identity
        of a set of hits between a query area and a reference area.

        Arguments:
            hits: an iterator of Hit instances
            count: the number of hits

        Returns:
            a float between 0 and 1
    """
    if not hits or not count:
        return 0.
    score = sum(hit.percent_identity for hit in hits) / count
    score /= 100  # from double digit percentage
    # avoid division by 0
    if score > 1. - 1e-8:
        return 1.
    # avoid a log of 0, though these *should* be culled already
    if score < 1e-8:
        return 0.
    # logit function, shifted from x range of (-6, 6) to (0, 1)
    # tiny values can result in small negatives
    # very high scores can result in just over 1 (e.g. .997 -> 1.003)
    result = min(1, max(0, math.log(score/(1-score)) / 12 + 0.5))
    assert 0 <= result <= 1
    return result


def trim_to_best_hit(hits_by_reference_gene: Dict[str, List[Hit]]) -> Dict[str, Hit]:
    """ Converts a N:M set of hits to a 1:1 set.
        Priority is based on percentage identity, any ties are resolved by ordering
        the reference CDS name.

        Arguments:
            hits_by_reference_gene: a dictionary mapping reference CDS name to a list of Hits

        Returns:
            a dictionary mapping reference CDS name to it's best hit
    """

    pairs = {}
    for ref, hits in hits_by_reference_gene.items():
        for hit in hits:
            pairs[(ref, hit)] = hit.percent_identity

    best_first = (i[0] for i in sorted(pairs.items(), key=lambda x: (x[1], x[0][0]), reverse=True))

    mapping: Dict[str, Hit] = {}

    features: Set[CDSFeature] = set()
    for ref, hit in best_first:
        if hit.cds in features or ref in mapping:
            continue
        mapping[ref] = hit
        assert hit.cds not in features
        features.add(hit.cds)
    assert 1 <= len(mapping) <= len(hits_by_reference_gene)
    return mapping


def convert_to_references(hits_by_name: HitsByReferenceName, references: Dict[str, ReferenceRecord]) -> HitsByReference:
    """ Converts results by reference accession and CDS numeric ID to results by
        ReferenceRegion and CDS name.

        Arguments:
            hits_by_name: a dictionary mapping reference accession to
                            a dictionary mapping reference CDS numeric ID to
                                a list of Hits for that reference
            references: a dictionary mapping reference accession to ReferenceRecord

        Returns:
            a dictionary mapping ReferenceRegion to
                a dictionary mapping reference CDS name to
                    a list of Hits for that reference
    """
    results: HitsByReference = defaultdict(dict)
    for name, hits in hits_by_name.items():
        reference_record = references[name]
        for cds_id, cds_hits in hits.items():
            cds_name = reference_record.cds_mapping[str(cds_id)]
            for region in reference_record.regions:
                assert region.accession == reference_record.accession
                if cds_name in region.cdses:
                    results[region][cds_name] = cds_hits
    return results


def find_diamond_matches(record: Record, database: str) -> Tuple[HitsByCDS, HitsByReferenceName]:
    """ Runs diamond, comparing all features in the record to the given database

        Arguments:
            record: the record to use as a query
            database: the path of the database to compare to

        Returns:
            a tuple of
                a dictionary mapping CDSFeature to
                     a dictionary mapping reference CDS numeric ID to
                        a list of Hits for that reference
                a dictionary mapping reference region name to
                    a dictionary mapping reference CDS numeric ID to
                        a list of Hits for that reference
    """
    logging.info("Comparing regions to reference database")
    extra_args = [
        "--compress", "0",
        "--max-target-seqs", "10000",
        "--evalue", "1e-05",
        "--outfmt", "6",  # 6 is blast tabular format, just as in blastp
    ]
    features = record.get_cds_features_within_regions()

    with NamedTemporaryFile() as temp_file:
        temp_file.write(fasta.get_fasta_from_features(features, numeric_names=True).encode())
        temp_file.flush()
        raw = subprocessing.run_diamond_search(temp_file.name, database, mode="blastp", opts=extra_args)
    return blast_parse(raw, dict(enumerate(features)))


def blast_parse(diamond_output: str, inputs_to_features: Dict[int, CDSFeature],
                min_seq_coverage: float = 30., min_perc_identity: float = 25.
                ) -> Tuple[HitsByCDS, HitsByReferenceName]:
    """ Parses diamond/blastp tabular output into a usable form.

        Arguments:
            diamond_output: the output from diamond in blast format
            inputs_to_features: a dictionary mapping numeric CDS index used in
                                the input FASTA to the matching CDSFeature
            min_seq_coverage: the exclusive lower bound of sequence coverage for a match
            min_perc_identity: the exclusive lower bound of identity similarity for a match

        Returns:
            a tuple of
                a dictionary mapping CDSFeature to
                     a dictionary mapping reference CDS numeric ID to
                        a list of Hits for that reference
                a dictionary mapping reference accession to
                    a dictionary mapping reference CDS numeric ID to
                        a list of Hits for that reference
    """
    hits_by_cds: HitsByCDS = defaultdict(lambda: defaultdict(list))
    hits_by_reference: HitsByReferenceName = defaultdict(lambda: defaultdict(list))
    for line in diamond_output.splitlines():
        hit = parse_hit(line, inputs_to_features)
        if not (hit.percent_identity >= min_perc_identity and hit.percent_coverage >= min_seq_coverage):
            continue
        hits_by_cds[hit.cds][hit.reference_record].append(hit)
        hits_by_reference[hit.reference_record][hit.reference_id].append(hit)
    return hits_by_cds, hits_by_reference


def parse_hit(line: str, feature_mapping: Dict[int, CDSFeature]) -> Hit:
    """ Parse a Hit instance from a diamond/blastp tabular output line

        Arguments:
            parts: the line to parse
            feature_mapping: a dictionary mapping numeric CDS index used in
                             the input FASTA to the matching CDSFeature

        Returns:
            a Hit instance
    """
    # blast tabular chunks
    # 0. 	 qseqid 	 query (e.g., gene) sequence id
    # 1. 	 sseqid 	 subject (e.g., reference genome) sequence id
    # 2. 	 pident 	 percentage of identical matches
    # 3. 	 length 	 alignment length
    # 4. 	 mismatch 	 number of mismatches
    # 5. 	 gapopen 	 number of gap openings
    # 6. 	 qstart 	 start of alignment in query
    # 7. 	 qend 	 end of alignment in query
    # 8. 	 sstart 	 start of alignment in subject
    # 9. 	 send 	 end of alignment in subject
    # 10. 	 evalue 	 expect value
    # 11. 	 bitscore 	 bit score

    parts = line.split("\t")
    assert len(parts) == 12

    cds = feature_mapping[int(parts[0])]
    reference_record, reference_id = parts[1].split("|")
    perc_ident = float(parts[2])
    evalue = float(parts[10])
    blastscore = float(parts[11])
    perc_coverage = (float(parts[3]) / len(cds.translation)) * 100

    return Hit(reference_record, reference_id, cds, perc_ident, blastscore, perc_coverage, evalue)
