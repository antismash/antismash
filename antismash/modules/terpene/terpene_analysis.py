# License: GNU Affero General Public License v3 or later
# A copy of GNU AGPL v3 should have been included in this software package in LICENSE.txt.

""" The analysis section of the terpene module
"""
import logging

from typing import Iterable
from collections import defaultdict

from antismash.common import fasta, path
from antismash.common.hmmscan_refinement import HMMResult, QueryResult, gather_by_query, remove_incomplete
from antismash.common.secmet import Record, Protocluster, CDSFeature
from antismash.common.subprocessing.hmmpfam import run_hmmpfam2
from antismash.common.subprocessing.hmmscan import run_hmmscan
from antismash.common.utils import get_query_positions_from_alignment

from .data_loader import (
    CompoundGroup, Motif, Reaction, Relationship, TerpeneHMM,
    load_hmm_lengths, load_hmm_properties, load_relationships,
    )
from .results import Activity, DomainPrediction, MotifResult, ProtoclusterPrediction


_MAIN_PROFILES_PATH =  path.get_full_path(__file__, "data", "main_profiles.hmm")


_MAIN_TYPE_PRIORITY = {val: i for i, val in enumerate([
    "Lycopene_cycl",
    "TS_Pyr4",
    "TS_UbiA",
    "T1TS",
    "HAD_2",
    "T2TS",
    "PT_phytoene_like",
    "PT_FPPS_like",
    "ambiguous"
])}


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


def group_hmm_results(hmm_results: list[HMMResult],
                      *, allowed_overlap: int = 20) -> list[list[HMMResult]]:
    """ Groups overlapping hmm results, if the overlap is larger than
        the allowed overlap

        Arguments:
            hmm_results: a list of HMMResult objects
            allowed_overlap: the maximum overlap of hits before creating a new group

        Returns:
            grouped_results: a list of lists of HMMResult objects

    """
    hmm_results = sorted(hmm_results, key=lambda result: result.query_start)
    groups = [[hmm_results[0]]]
    end = hmm_results[0].query_end
    for hmm_result in hmm_results[1:]:
        if hmm_result.query_start < (end - allowed_overlap):
            groups[-1].append(hmm_result)
        else:
            groups.append([hmm_result])
        end = max(hmm_result.query_end, end)
    return groups


def filter_subtypes(hmm_results: list[HMMResult],
                    relationships: dict[str, Relationship],
                    check_parents: bool = False) -> list[HMMResult]:
    """ Removes hits if a more specific child hit is present

        Arguments:
            hmm_results: a list of HMMResult objects
            relationships: a dictionary of hmm names to Relationship objects
            check_parents: if True, removes hits when a hit to the parent profile is absent

        Returns:
            a list of HMMResult objects
    """
    remaining = []
    results_by_name = {result.hit_id : result for result in hmm_results}
    chains_by_name : dict[str, list[list]] = {}

    def get_ancestor_chains(hmm_name: str) -> list[list]:
        parents = relationships[hmm_name].parents
        if not parents:
            return [[hmm_name]]
        chains = []
        for parent in parents:
            for chain in get_ancestor_chains(parent):
                chains.append([hmm_name] + chain)
        return chains

    for result in hmm_results:
        if check_parents:
            try:
                chains = chains_by_name[result.hit_id]
            except KeyError:
                full_chains = get_ancestor_chains(result.hit_id)
                chains = [chain[1:] for chain in full_chains]
                chains_by_name[result.hit_id] = chains
                valid = False
                for chain in chains:
                    if all(ancestor in results_by_name
                            and result.overlaps_with(results_by_name[ancestor])
                            for ancestor in chain):
                        valid = True
                if not valid:
                    logging.debug("%s: Missing higher-level hit, dropping result.", result)
                    continue
        children = relationships[result.hit_id].children
        if any(child in results_by_name for child in children):
            continue
        remaining.append(result)
    return remaining


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


def extract_motif(query: str, hit: QueryResult, offset: int,
                  motif: Motif) -> MotifResult:
    """ Extracts a motif from a query sequence based on the alignment to a profile, and evaluates
        its activity based on whether the extracted residues match the expected residues for the motif.

        Arguments:
            query: the query sequence
            hit: the alignment of the query to the profile
            offset: offset of the query compared to the original sequence (for calculating absolute positions)
            motif: the motif to extract

        Returns:
            a MotifResult object
    """
    query_positions = get_query_positions_from_alignment(hit, motif.positions)
    if not query_positions:
        return MotifResult(name=motif.name, active=False)
    query_residues = [query[i] for i in query_positions]
    # Check if the positions match the expected amino acids for the motif
    active = True
    for query_residue, allowed_residues in zip(query_residues, motif.allowed_residues):
        if not query_residue in allowed_residues:
            active = False
            break
    start = query_positions[0]
    end = query_positions[-1]+1 # Python-style exclusive end position
    motif_seq = query[start:end]
    # Start and end positions of the motif are stored relative to the original sequence
    motif_start = start + offset
    motif_end = end + offset
    return MotifResult(name=motif.name, active=active, sequence=motif_seq, start=motif_start, end=motif_end)


def get_motif_results(hmm_file: str, domain_seq: str, domain_start: int,
                      hmm_properties: dict[str, TerpeneHMM], domain_type: str,
                      ) -> tuple[MotifResult, ...]:
    """ Extracts any active site motifs present in a domain sequence and evaluates their activity.

        Arguments:
            hmm_file: the path to an HMM file to use for motif extraction
            domain_seq: the sequence of the domain
            domain_start: the start position of the domain within the full CDS sequence
            hmm_properties: a dictionary of hmm names to TerpeneHMM objects
            domain_type: the type of the domain, used to select the correct profile for motif extraction

        Returns:
            a tuple of MotifResult objects
    """
    hmmpfam_results = run_hmmpfam2(hmm_file, f">query\n{domain_seq}\n")
    hits = hmmpfam_results[0].hsps
    for hit in sorted(hits, key=lambda h: h.bitscore, reverse=True):
        hit_type = hit.hit.id
        if domain_type in hit_type: # Covers both T1TS_KS and T1TS profiles for T1TS domains
            profile = hmm_properties[hit_type]
            motif_results = []
            for motif in profile.motifs:
                motif_results.append(extract_motif(domain_seq, hit, domain_start, motif))
            return tuple(motif_results)
    return tuple()


def get_domain_activity(motif_results: tuple[MotifResult, ...]) -> Activity:
    """ Determines the overall activity of a domain based on its motif results.
        All motifs must be active for the domain to be considered active.

        Arguments:
            motif_results: a tuple of MotifResult objects

        Returns:
            an Activity enum value
    """
    if not motif_results:
        return Activity.UNKNOWN
    if all(result.active for result in motif_results):
        return Activity.ACTIVE
    return Activity.INACTIVE


def get_domain_prediction(hmm_results: list[HMMResult],
                          hmm_properties: dict[str, TerpeneHMM],
                          relationships: dict[str, Relationship],
                          cds_translation: str) -> DomainPrediction:
    """ Converts a list of HMMResults to a DomainPrediction

        Arguments:
            hmm_results: a list of HMMResults

        Returns:
            domain_prediction: a DomainPrediction
    """
    start = min(hmm_result.query_start for hmm_result in hmm_results)
    end = max(hmm_result.query_end for hmm_result in hmm_results)
    hmm_results = filter_subtypes(hmm_results, relationships)
    profiles = [hmm_properties[hmm_result.hit_id] for hmm_result in hmm_results]
    main_types = set(profile.domain_type for profile in profiles)
    subtypes = set(profile.name for profile in profiles if profile.is_subtype())
    if not hmm_results or len(main_types) > 1:
        logging.debug("Ambiguous hit detected.")
        return DomainPrediction(domain_type="ambiguous", subtypes=tuple(),
                                start=start, end=end,
                                reactions=tuple(),
                                motif_results=tuple(), activity=Activity.UNKNOWN)

    domain_type = main_types.pop()
    final_reactions = merge_reactions_by_substrate(profiles)

    # Extract active site motifs and their activity
    hmm_file = _MAIN_PROFILES_PATH
    domain_seq = cds_translation[start:end]
    assert len(domain_seq) == end-start, "Domain indexes out of range"

    motif_results = get_motif_results(hmm_file, domain_seq, start, hmm_properties, domain_type)

    # Get overall activity of the domain
    activity = get_domain_activity(motif_results)

    return DomainPrediction(domain_type=domain_type, subtypes=tuple(subtypes),
                            start=start, end=end,
                            reactions=final_reactions,
                            motif_results=motif_results, activity=activity)


def get_cds_predictions(hmm_results_per_cds: dict[str, list[HMMResult]],
                        hmm_properties: dict[str, TerpeneHMM],
                        relationships: dict[str, Relationship],
                        cds_by_name: dict[str, CDSFeature]) -> dict[str, list[DomainPrediction]]:
    """ Convert list of HMMResults in CDS mapping to a list of DomainPredictions

        Arguments:
            hmm_results: a mapping of CDS name to a list of HMMResults

        Returns:
            cds_predictions: a mapping of CDS name to a list of DomainPredictions
    """
    cds_predictions: dict[str, list[DomainPrediction]] = defaultdict(list)

    for cds_name, hmm_results in hmm_results_per_cds.items():
        cds_translation = cds_by_name[cds_name].translation
        grouped_results = group_hmm_results(hmm_results)
        # Add a domain prediction for each group
        for group in grouped_results:
            domain = get_domain_prediction(group, hmm_properties, relationships, cds_translation)
            cds_predictions[cds_name].append(domain)
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
    # order by whether the domain has any reactions
    # then the predefined main type priority
    ordered_domains = sorted(all_domains, key=lambda domain: (
        0 if domain.reactions else 1,
        _MAIN_TYPE_PRIORITY[domain.domain_type],
    ))

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


def analyse_cluster(cluster: Protocluster, record: Record) -> ProtoclusterPrediction:
    """ Analyse a terpene cluster

        Arguments:
            cluster: the Protocluster to analyse

        Returns:
            a single ProtoclusterPrediction instance with analysis results
    """
    assert cluster.product_category == "terpene"
    hmm_properties = load_hmm_properties()
    hmm_lengths = load_hmm_lengths(hmm_properties)
    relationships = load_relationships(hmm_properties)

    hmmscan_results = run_terpene_hmmscan(cluster.cds_children)
    refined_results = filter_incomplete(hmmscan_results, hmm_lengths)
    refined_results = filter_by_score(refined_results, hmm_properties)
    cds_predictions = get_cds_predictions(refined_results, hmm_properties,
                                          relationships, record.get_cds_name_mapping())
    return get_cluster_prediction(cds_predictions)
