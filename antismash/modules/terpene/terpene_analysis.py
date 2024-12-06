# License: GNU Affero General Public License v3 or later
# A copy of GNU AGPL v3 should have been included in this software package in LICENSE.txt.

""" The analysis section of the terpene module
"""
import logging

from typing import Iterable
from collections import defaultdict

from antismash.common import fasta, path
from antismash.common.hmmscan_refinement import HMMResult, QueryResult, gather_by_query, remove_incomplete
from antismash.common.secmet import Protocluster, CDSFeature
from antismash.common.subprocessing.hmmscan import run_hmmscan

from .data_loader import CompoundGroup, Reaction, TerpeneHMM, load_hmm_lengths, load_hmm_properties
from .results import ProtoclusterPrediction, DomainPrediction


_MAIN_TYPES_PRIORITY = (
    "Lycopene_cycl",
    "TS_Pyr4",
    "TS_UbiA",
    "T1TS",
    "HAD_2",
    "T2TS",
    "PT_phytoene_like",
    "PT_FPPS_like",
    "ambiguous"
)


def run_terpene_hmmscan(cds_features: Iterable[CDSFeature]) -> list[QueryResult]:
    """ Runs hmmscan for terpene proteins on the given CDSFeatures

        Arguments:
            cluster: Protocluster on which to run the terpene hmmscan

        Returns:
            hmmscan_results: a list of QueryResult objects from Bio's SearchIO
    """
    cluster_fasta = fasta.get_fasta_from_features(cds_features)
    hmm_file = path.get_full_path(__file__, "data", "all_profiles.hmm")
    hmmscan_results = run_hmmscan(target_hmmfile=hmm_file, query_sequence=cluster_fasta)

    return hmmscan_results


def filter_incomplete(hmmscan_results: list[QueryResult], hmm_lengths: dict[str, int]
                      ) -> dict[str, list[HMMResult]]:
    """ Removes incomplete fragments

        Arguments:
            hmmscan_results: a list of QueryResult objects from Bio's SearchIO
            hmm_lengths: a dictionary of hmm names to hmm lengths

        Returns:
            refined_results: a dictionary of cds names to HMMResult objects
    """
    results_by_id = gather_by_query(hmmscan_results)
    refined_results = {}
    for cds, results in results_by_id.items():
        refined = sorted(list(results), key=lambda result: result.query_start)
        refined = remove_incomplete(refined, hmm_lengths)
        if refined:
            refined_results[cds] = refined
    return refined_results


def filter_by_score(hmm_results_per_cds: dict[str, list[HMMResult]],
                    hmm_properties: dict[str, TerpeneHMM]) -> dict[str, list[HMMResult]]:
    """ Removes hmm results with bitscores that are below the specified cutoff for the hmm

        Arguments:
            hmm_results_per_cds: a dictionary of cds names to HMMResult objects
            hmm_properties: a dictionary of hmm names to TerpeneHMM objects

        Returns:
            refined_results: a dictionary of cds names to HMMResult objects

    """
    results_by_id: dict[str, list[HMMResult]] = defaultdict(list)
    for cds_name, results in hmm_results_per_cds.items():
        # Store result if it is above the cut-off
        for result in results:
            if result.hit_id in hmm_properties:
                terpene_hmm = hmm_properties[result.hit_id]
            else:
                raise ValueError(f"Failed to find signature for ID {result.hit_id}")
            if result.bitscore > terpene_hmm.cutoff:
                results_by_id[cds_name].append(result)
    return results_by_id


def group_hmm_results(hmm_results: list[HMMResult]) -> list[list[HMMResult]]:
    """ Groups overlapping hmm results, if the overlap is larger than
        the allowed overlap

        Arguments:
            hmm_results: a list of HMMResult objects

        Returns:
            grouped_results: a list of lists of HMMResult objects

    """
    hmm_results = sorted(hmm_results, key=lambda result: result.query_start)
    allowed_overlap = 20
    groups = [[hmm_results[0]]]
    end = hmm_results[0].query_end
    for hmm_result in hmm_results[1:]:
        if hmm_result.query_start < (end - allowed_overlap):
            groups[-1].append(hmm_result)
        else:
            groups.append([hmm_result])
        end = max(hmm_result.query_end, end)
    return groups


def merge_reactions_by_substrate(profiles: list[TerpeneHMM]
                                 ) -> tuple[Reaction, ...]:
    """ Merges the reactions for a group of hmms.
        Only merges reactions with identical substrates in
        all members of the group.
        If no reactions fit this criterium, outputs an empty tuple.

        Arguments:
            hmm_results: a list of HMMResult objects

        Returns:
            reactions: a tuple of Reaction objects

    """
    if not profiles:
        return tuple()
    # gather all reactions from all groups with the same substrate
    reactions_by_substrate: dict[tuple[CompoundGroup, ...], list[Reaction]] = defaultdict(list)
    profiles_with_reactions: int = 0
    for profile in profiles:
        if profile.reactions:
            for reaction in profile.reactions:
                reactions_by_substrate[reaction.substrates].append(reaction)
            profiles_with_reactions += 1
    # find the substrates which are present in all groups
    intersecting_groups = {sub: reactions for sub, reactions in reactions_by_substrate.items()
                           if len(reactions) == profiles_with_reactions}
    # if no substrates are present in all groups, return nothing
    if not intersecting_groups:
        return tuple()

    # for each matching substrate, compress the reaction down to the intersection of all products
    results = []
    for substrate_group in intersecting_groups.values():
        first = substrate_group[0]
        for reaction in substrate_group[1:]:
            first = first.build_intersection(reaction)
        # if no products remain, don't include the substrate
        if first.products:
            results.append(first)

    # return all remaining merged reactions
    return tuple(results)


def get_domain_prediction(hmm_results: list[HMMResult],
                          hmm_properties: dict[str, TerpeneHMM]) -> DomainPrediction:
    """ Converts a list of HMMResults to a DomainPrediction

        Arguments:
            hmm_results: a list of HMMResults

        Returns:
            domain_prediction: a DomainPrediction
    """
    start = min(hmm_result.query_start for hmm_result in hmm_results)
    end = max(hmm_result.query_end for hmm_result in hmm_results)
    profiles = [hmm_properties[hmm_result.hit_id] for hmm_result in hmm_results]
    main_types = set(profile.domain_type for profile in profiles)
    subtypes = set(profile.name for profile in profiles if profile.is_subtype())
    if len(main_types) > 1:
        logging.debug("Overlapping hits for different main types.")
        domain_type = "ambiguous"
        subtypes = set()
    else:
        domain_type = main_types.pop()
    final_reactions = merge_reactions_by_substrate(profiles)
    return DomainPrediction(domain_type=domain_type, subtypes=tuple(subtypes),
                            start=start, end=end,
                            reactions=final_reactions)


def get_cds_predictions(hmm_results_per_cds: dict[str, list[HMMResult]],
                        hmm_properties: dict[str, TerpeneHMM]) -> dict[str, list[DomainPrediction]]:
    """ Convert list of HMMResults in CDS mapping to a list of DomainPredictions

        Arguments:
            hmm_results: a mapping of CDS name to a list of HMMResults

        Returns:
            cds_predictions: a mapping of CDS name to a list of DomainPredictions
    """
    cds_predictions: dict[str, list[DomainPrediction]] = defaultdict(list)

    for cds_name, hmm_results in hmm_results_per_cds.items():
        grouped_results = group_hmm_results(hmm_results)
        # Add a domain prediction for each group
        for group in grouped_results:
            cds_predictions[cds_name].append(get_domain_prediction(group, hmm_properties))
    return cds_predictions


def get_cluster_prediction(cds_predictions: dict[str, list[DomainPrediction]]) -> ProtoclusterPrediction:
    """ Convert list of DomainPredictions in CDS mapping to a ProtoclusterPrediction

        Arguments:
            cds_predictions: a mapping of CDS name to a list of DomainPredictions

        Returns:
            a single ProtoclusterPrediction instance
    """
    cluster_pred = ProtoclusterPrediction(cds_predictions)
    all_domains = [domain for domains in cds_predictions.values() for domain in domains]
    if not all_domains:
        return cluster_pred
    # Order domains by the predefined main type priority
    ordered_domains = sorted(all_domains, key=lambda domain: _MAIN_TYPES_PRIORITY.index(domain.domain_type))
    # Then order by whether the domain has any reactions
    ordered_domains = sorted(ordered_domains, key=lambda domain: not domain.reactions)

    # Find the products of the most informative domain(s)
    compounds_by_name: dict[str, set[CompoundGroup]] = defaultdict(set)
    domain_type = ordered_domains[0].domain_type
    for domain in ordered_domains:
        if domain.domain_type != domain_type:
            break
        for reaction in domain.reactions:
            for product in reaction.products:
                compounds_by_name[product.name].add(product)
    for compounds in sorted(compounds_by_name.values()):
        assert len(compounds) == 1, "Multiple CompoundGroup instances with the same name are not allowed"
        cluster_pred.products.append(compounds.pop())
    return cluster_pred


def analyse_cluster(cluster: Protocluster) -> ProtoclusterPrediction:
    """ Analyse a terpene cluster

        Arguments:
            cluster: the Protocluster to analyse

        Returns:
            a single ProtoclusterPrediction instance with analysis results
    """
    assert cluster.product_category == "terpene"
    hmm_properties = load_hmm_properties()
    hmm_lengths = load_hmm_lengths(hmm_properties)

    hmmscan_results = run_terpene_hmmscan(cluster.cds_children)
    refined_results = filter_incomplete(hmmscan_results, hmm_lengths)
    refined_results = filter_by_score(refined_results, hmm_properties)
    cds_predictions = get_cds_predictions(refined_results, hmm_properties)
    return get_cluster_prediction(cds_predictions)
