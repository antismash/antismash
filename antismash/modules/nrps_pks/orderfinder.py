# License: GNU Affero General Public License v3 or later
# A copy of GNU AGPL v3 should have been included in this software package in LICENSE.txt.

""" Calculates a likely order of NRPS/PKS domains """

import itertools
import logging
import os
from typing import Dict, List, Optional, Tuple

from antismash.common import path, subprocessing, utils
from antismash.common.secmet import CDSFeature, Record

from .results import SuperClusterPrediction


def analyse_biosynthetic_order(nrps_pks_features: List[CDSFeature],
                               consensus_predictions: Dict[str, str],
                               record: Record) -> List[SuperClusterPrediction]:
    """ For each NRPS or PKS supercluster, determines if that supercluster is
        docking or not then determines the monomer ordering

        Arguments:
            nrps_pks_features: all NRPS/PKS features within the record
            consensus_predictions: a dictionary mapping each NRPS/PKS domain name to its prediction
            record: the Record being analysed

        Returns:
            a dictionary mapping supercluster number to
                a tuple of
                    prediction string
                    and whether docking domain analysis was used for the prediction
    """
    compound_predictions = []  # type: List[SuperClusterPrediction]
    # Find NRPS/PKS gene superclusters
    superclusters = [cluster for cluster in record.get_superclusters()
                             if "nrps" in "-".join(cluster.products)
                             or "pks" in "-".join(cluster.products)]
    if not superclusters:
        return []
    # Predict biosynthetic gene order in superclusters using starter domains,
    # thioesterase domains, gene order and docking domains
    for supercluster in superclusters:
        supercluster_number = supercluster.get_supercluster_number()
        cds_in_supercluster = [gene for gene in nrps_pks_features if gene.overlaps_with(supercluster)]
        if not cds_in_supercluster:
            continue
        pks_features, nrps_count, hybrid_count = find_supercluster_modular_enzymes(cds_in_supercluster)
        # If more than three PKS cds features, use dock_dom_analysis if possible to identify order
        # since this will grow as n!, an upper limit is also required
        if 3 < len(pks_features) < 11 and not nrps_count and not hybrid_count:
            logging.debug("SuperCluster %d monomer ordering method: domain docking analysis", supercluster_number)
            geneorder = perform_docking_domain_analysis(pks_features)
            docking = True
        else:
            logging.debug("SuperCluster %d monomer ordering method: colinear", supercluster_number)
            geneorder = find_colinear_order(cds_in_supercluster)
            docking = False

        prediction = generate_substrates_order(geneorder, consensus_predictions)
        compound_predictions.append(SuperClusterPrediction(supercluster_number, prediction, docking))
    return compound_predictions


def find_supercluster_modular_enzymes(cds_features: List[CDSFeature]) -> Tuple[List[CDSFeature], int, int]:
    """ Finds PKS-only features and counts the NRPS-only features and
        hybrid (not a combination of individual PKS and NRPS domains)
        features in a set of CDS features.

        Arguments:
            cds_features: the CDS features to process

        Returns:
            a tuple of
                a list of PKS-only CDS features
                a count of NRPS-only features
                a count of hybrid features
    """
    pks_features = []
    nrps_count = 0
    hybrid_count = 0
    for cds in cds_features:
        classification = cds.nrps_pks.type
        if "PKS" in classification and "NRPS" not in classification:
            pks_features.append(cds)
        elif "PKS" not in classification and "NRPS" in classification:
            nrps_count += 1
        elif "PKS/NRPS" in classification:
            domain_names = set(cds.nrps_pks.domain_names)
            contains_nrps = domain_names.intersection({"AMP-binding", "A-OX", "Condensation"})
            contains_pks = domain_names.intersection({"PKS_KS", "PKS_AT"})
            if contains_pks and not contains_nrps:
                pks_features.append(cds)
            # the case of both single pks domain and nrps domain(s) is ignored
            # because that construction isn't meaningful
        elif "Hybrid" in classification:
            hybrid_count += 1
    return pks_features, nrps_count, hybrid_count


def generate_substrates_order(geneorder: List[CDSFeature], consensus_predictions: Dict[str, str]) -> str:
    """ Generate substrates order from predicted gene order and consensus
        predictions. E.g. (ala-dpg) + (pk).

        Arguments:
            geneorder: a list of CDSFeatures
            consensus_predictions: a dictionary mapping domain name to prediction

        Returns:
            the resulting monomer prediction string
    """
    predictions = []

    for gene in geneorder:
        consensuses = []
        for domain in gene.nrps_pks.domains:
            consensus = consensus_predictions.get(domain.feature_name)
            if consensus:
                consensuses.append(consensus)
        if consensuses:
            predictions.append("(%s)" % (" - ".join(consensuses)))

    if not predictions:
        return ""
    return " + ".join(predictions)


def find_first_and_last_cds(cds_features: List[CDSFeature]) -> Tuple[Optional[CDSFeature], Optional[CDSFeature]]:
    """ Find first and last CDSFeature based on starter module and TE / TD.

        If multiple possibilities are found for start or end, no gene will be
        returned as such.

        Arguments:
            cds_features: the CDS features to search in for start and end CDS

        Returns:
            a tuple of
                the start CDS or None, and
                the end CDS or None
    """

    start_cds = None
    end_cds = None

    # find the end
    for cds in cds_features:
        domain_names = cds.nrps_pks.domain_names
        if "Thioesterase" in domain_names or "TD" in domain_names:
            if end_cds:
                end_cds = None
                break
            end_cds = cds

    # find the start
    for cds in cds_features:
        if cds == end_cds:
            continue
        domain_names = cds.nrps_pks.domain_names
        if domain_names[:2] == ["PKS_AT", "ACP"]:
            if start_cds:
                # two possible starts, don't attempt fallbacks
                return None, end_cds
            start_cds = cds

    # if no AT-ACP start gene, try looking for KS-AT-ACP
    if not start_cds:
        for cds in cds_features:
            if cds == end_cds:
                continue
            domain_names = cds.nrps_pks.domain_names
            if domain_names[:3] == ["PKS_KS", "PKS_AT", "ACP"]:
                if start_cds:
                    start_cds = None
                    break
                start_cds = cds
    return start_cds, end_cds


def extract_nterminus(data_dir: str, cds_features: List[CDSFeature], start_cds: Optional[CDSFeature]) -> Dict[str, str]:
    """ -extract N-terminal 50 residues of each non-starting protein
        -scan for docking domains using hmmsearch
        -parse output to locate interacting residues
    """
    n_terminal_residues = {}
    n_terminals = {}
    nterm_file = os.path.join(data_dir, 'nterm.fasta')
    for cds in cds_features:
        if cds is not start_cds:
            seq = str(cds.translation)
            n_terminals[cds.get_name()] = seq[:50]
    for name, seq in n_terminals.items():
        alignments = subprocessing.run_muscle_single(name, seq, nterm_file)
        query_seq = alignments[name]
        ref_seq = alignments["EryAIII_5_6_ref"]
        n_terminal_residues[name] = utils.extract_by_reference_positions(query_seq, ref_seq, [2, 15])
    return n_terminal_residues


def extract_cterminus(data_dir: str, cds_features: List[CDSFeature], end_cds: Optional[CDSFeature]) -> Dict[str, str]:
    """ Extract C-terminal 100 residues of each non-ending protein,
        scan for docking domains, parse output to locate interacting residues

        Arguments:
            data_dir: the directory containing the C-terminal reference files
            cds_features: the list of CDSFeatures to extract terminals from
            end_cds: if not None, skips this CDS since C-terminals are irrelevant

        Returns:
            A dictionary mapping gene name to the pair of residues extracted
    """
    c_terminal_residues = {}
    c_terminals = {}  # type: Dict[str, str]
    cterm_file = os.path.join(data_dir, 'cterm.fasta')
    for cds in cds_features:
        if cds is not end_cds:
            seq = str(cds.translation)
            c_terminals[cds.get_name()] = seq[-100:]
    for name, seq in c_terminals.items():
        alignments = subprocessing.run_muscle_single(name, seq, cterm_file)
        query_seq = alignments[name]
        ref_seq = alignments["EryAII_ref"]
        c_terminal_residues[name] = utils.extract_by_reference_positions(query_seq, ref_seq, [55, 64])
    return c_terminal_residues


def find_possible_orders(cds_features: List[CDSFeature], start_cds: Optional[CDSFeature],
                         end_cds: Optional[CDSFeature]) -> List[List[CDSFeature]]:
    """ Finds all possible arrangements of the given cds_features. If not None, the
        start gene will always be the first in each order. Similarly, the end
        gene will always be last.

        Arguments:
            cds_features: a list of all CDSFeatures, may include start_cds and end_cds
            start_cds: None or the CDS with which to start every arrangement
            end_cds: None or the CDS with which to end every arrangement

        Returns:
            a list of lists, each sublist being a unique ordering of the
            provided CDSFeatures
    """
    assert len(cds_features) < 11, "input too large, function is O(n!)"
    assert start_cds is None or isinstance(start_cds, CDSFeature)
    assert end_cds is None or isinstance(end_cds, CDSFeature)
    if start_cds or end_cds:
        assert start_cds != end_cds, "Using same gene for start and end of ordering"
    cds_to_order = []
    for cds in cds_features:
        if cds == start_cds or cds == end_cds:
            pass
        else:
            cds_to_order.append(cds)
    possible_orders = []
    start = []  # type: List[CDSFeature]
    if start_cds:
        start = [start_cds]
    end = []  # type: List[CDSFeature]
    if end_cds:
        end = [end_cds]
    for order in itertools.permutations(cds_to_order, len(cds_to_order)):
        possible_orders.append(start + list(order) + end)
    # ensure the list of possible orders is itself ordered for reliability
    return sorted(possible_orders, key=lambda x: [g.location.start for g in x])


def rank_biosynthetic_orders(n_terminal_residues: Dict[str, str],
                             c_terminal_residues: Dict[str, str],
                             possible_orders: List[List[CDSFeature]]) -> List[CDSFeature]:
    """ Scores each possible order according to terminal pairs of adjacent cds_features.

        Arguments:
            n_terminal_residues: a dictionary mapping CDSFeature to their pair of N terminal residues
            c_terminal_residues: a dictionary mapping CDSFeature to their pair of C terminal residues
            possible_orders: a list of gene orderings to evaluate

        Returns:
            the first ordering that scored highest or equal highest
    """
    assert possible_orders
    # If docking domains found in all, check for optimal order using interacting residues
    hydrophobic = {"A", "V", "I", "L", "F", "W", "Y", "M"}
    positively_charged = {"H", "K", "R"}
    negatively_charged = {"D", "E"}
    # find best scoring order
    best_score = -2 * len(possible_orders[0])
    best_order = possible_orders[0]
    for order in possible_orders:
        score = 0
        interactions = [order[i:i + 2] for i in range(len(order) - 1)]
        for gene, next_gene in interactions:
            res1a, res2a = c_terminal_residues[gene.get_name()]
            res1b, res2b = n_terminal_residues[next_gene.get_name()]
            for pair in [{res1a, res1b}, {res2a, res2b}]:
                both_hydrophobic = pair.issubset(hydrophobic)
                same_polarity = pair.issubset(positively_charged) or pair.issubset(negatively_charged)
                opposite_polarity = len(pair & positively_charged) * len(pair & negatively_charged) == 1
                if both_hydrophobic or opposite_polarity:
                    score += 1
                elif same_polarity:
                    score -= 1
        if score > best_score:
            best_order = order
            best_score = score
    return best_order


def perform_docking_domain_analysis(cds_features: List[CDSFeature]) -> List[CDSFeature]:
    """ Estimates gene ordering based on docking domains of features

        Arguments:
            cds_features: a list of CDSFeatures to order

        Returns:
            a list of CDSFeatures in estimated order
    """
    start_cds, end_cds = find_first_and_last_cds(cds_features)
    data_dir = path.get_full_path(__file__, "data", "terminals")

    n_terminal_residues = extract_nterminus(data_dir, cds_features, start_cds)
    c_terminal_residues = extract_cterminus(data_dir, cds_features, end_cds)
    possible_orders = find_possible_orders(cds_features, start_cds, end_cds)

    geneorder = rank_biosynthetic_orders(n_terminal_residues, c_terminal_residues, possible_orders)
    return geneorder


def find_colinear_order(cds_features: List[CDSFeature]) -> List[CDSFeature]:
    """ Estimates gene ordering based on colinearity

        Arguments:
            cds_features: a list of CDSFeatures to order

        Returns:
            a list of CDSFeatures in estimated order
    """
    direction = 0
    for gene in cds_features:
        direction += gene.strand
    geneorder = list(cds_features)
    # Reverse if first gene encodes a multidomain protein with a TE/TD domain
    if direction < 0:
        geneorder.reverse()
    gene_domains = geneorder[0].nrps_pks.domain_names
    if "Thioesterase" in gene_domains or "TD" in gene_domains:
        if len(gene_domains) > 1:
            geneorder.reverse()
    return geneorder
