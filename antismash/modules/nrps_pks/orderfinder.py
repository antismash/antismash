# License: GNU Affero General Public License v3 or later
# A copy of GNU AGPL v3 should have been included in this software package in LICENSE.txt.

""" Calculates a likely order of NRPS/PKS domains """

import itertools
import logging
import os
from typing import Dict, List, Optional, Tuple

from antismash.common import path, subprocessing, utils
from antismash.common.secmet import CDSFeature


def analyse_biosynthetic_order(nrps_pks_genes, consensus_predictions, seq_record) -> Dict[int, Tuple[str, bool]]:
    """ For each NRPS or PKS cluster, determines if that cluster is docking or not
        then calls generate_substrates_order()

        Arguments:
            gene_domains: a dict mapping gene name to list of domains in that gene

        Returns:
            a dictionary mapping cluster number to
                a tuple of
                    prediction string
                    and whether docking domain analysis was used for the prediction
    """
    compound_predictions = {}  # type: Dict[int, Tuple[str, bool]]
    # Find NRPS/PKS gene clusters
    nrpspksclusters = [cluster for cluster in seq_record.get_clusters()
                       if "nrps" in cluster.products or "pks" in "-".join(cluster.products)]
    if not nrpspksclusters:
        return {}
    # Predict biosynthetic gene order in gene cluster using starter domains,
    # thioesterase domains, gene order and docking domains
    for cluster in nrpspksclusters:
        cluster_number = cluster.get_cluster_number()
        genes_in_cluster = [gene for gene in nrps_pks_genes if gene.overlaps_with(cluster)]
        if not genes_in_cluster:
            continue
        pks_count, nrps_count, hybrid_count = find_cluster_modular_enzymes(genes_in_cluster)
        # If more than three PKS genes, use dock_dom_analysis if possible to identify order
        if 3 < pks_count < 11 and not nrps_count and not hybrid_count:
            logging.debug("Cluster %d monomer ordering method: domain docking analysis", cluster_number)
            geneorder = perform_docking_domain_analysis(genes_in_cluster)
            docking = True
        else:
            logging.debug("Cluster %d monomer ordering method: colinear", cluster_number)
            geneorder = find_colinear_order(genes_in_cluster)
            docking = False
        prediction = generate_substrates_order(geneorder, consensus_predictions)
        compound_predictions[cluster_number] = (prediction, docking)
    return compound_predictions


def find_cluster_modular_enzymes(genes) -> Tuple[int, int, int]:
    """ counts number of PKS domains, NRPS domains and hybrid domains in a cluster
    """
    pksgenes = 0
    nrpsgenes = 0
    hybridgenes = 0
    for gene in genes:
        classification = gene.nrps_pks.type
        if "PKS" in classification and "NRPS" not in classification:
            pksgenes += 1
        elif "PKS" not in classification and "NRPS" in classification:
            nrpsgenes += 1
        elif "PKS/NRPS" in classification:
            domain_names = set(gene.nrps_pks.domain_names)
            contains_nrps = domain_names.intersection({"AMP-binding", "A-OX", "Condensation"})
            contains_pks = domain_names.intersection({"PKS_KS", "PKS_AT"})
            if contains_pks and not contains_nrps:
                pksgenes += 1
            # the case of both single pks domain and nrps domain(s) is ignored
            # because that construction isn't meaningful
        elif "Hybrid" in classification:
            hybridgenes += 1
    return pksgenes, nrpsgenes, hybridgenes


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
        gene_name = gene.get_name()
        for domain in gene.nrps_pks.domains:
            domain_name = gene_name + domain.label
            consensus = consensus_predictions.get(domain_name)
            if consensus:
                consensuses.append(consensus)
        if consensuses:
            predictions.append("(%s)" % ("-".join(consensuses)))

    if not predictions:
        return ""
    return " + ".join(predictions)


def find_first_and_last_genes(genes: List[CDSFeature]) -> Tuple[Optional[CDSFeature], Optional[CDSFeature]]:
    """ Find first and last genes based on starter module and TE / TD.

        If multiple possibilities are found for start or end, no gene will be
        returned as such.

        Arguments:
            genes: the CDS features to search in for start and end genes

        Returns:
            a tuple of
                the start gene or None, and
                the end gene or None
    """

    start_gene = None
    end_gene = None

    # find the end
    for gene in genes:
        domain_names = gene.nrps_pks.domain_names
        if "Thioesterase" in domain_names or "TD" in domain_names:
            if end_gene:
                end_gene = None
                break
            end_gene = gene

    # find the start
    for gene in genes:
        domain_names = gene.nrps_pks.domain_names
        if domain_names[:2] == ["PKS_AT", "ACP"]:
            if start_gene:
                # two possible starts, don't attempt fallbacks
                return None, end_gene
            start_gene = gene

    # if no AT-ACP start gene, try looking for KS-AT-ACP
    if not start_gene:
        for gene in genes:
            domain_names = gene.nrps_pks.domain_names
            if domain_names[:3] == ["PKS_KS", "PKS_AT", "ACP"]:
                if start_gene:
                    start_gene = None
                    break
                start_gene = gene
    return start_gene, end_gene


def extract_nterminus(data_dir, genes, start_gene):
    """ -extract N-terminal 50 residues of each non-starting protein
        -scan for docking domains using hmmsearch
        -parse output to locate interacting residues
    """
    n_terminal_residues = {}
    n_terminals = {}
    nterm_file = os.path.join(data_dir, 'nterm.fasta')
    for gene in genes:
        gene_name = gene.get_name()
        if gene_name != start_gene:
            seq = str(gene.translation)
            n_terminals[gene_name] = seq[:50]
    for name, seq in n_terminals.items():
        alignments = subprocessing.run_muscle_single(name, seq, nterm_file)
        query_seq = alignments[name]
        ref_seq = alignments["EryAIII_5_6_ref"]
        n_terminal_residues[name] = utils.extract_by_reference_positions(query_seq, ref_seq, [2, 15])
    return n_terminal_residues


def extract_cterminus(data_dir, genes, end_gene) -> Dict[str, str]:
    """ Extract C-terminal 100 residues of each non-ending protein,
        scan for docking domains, parse output to locate interacting residues

        Arguments:
            data_dir: the directory containing the C-terminal reference files
            genes: the list of genes to extract terminals from
            end_gene: if not None, skips this gene since C-terminals are irrelevant

        Returns:
            A dictionary mapping gene name to the pair of residues extracted
    """
    c_terminal_residues = {}
    c_terminals = {}  # type: Dict[str, str]
    cterm_file = os.path.join(data_dir, 'cterm.fasta')
    for gene in genes:
        gene_name = gene.get_name()
        if gene_name != end_gene:
            seq = str(gene.translation)
            c_terminals[gene_name] = seq[-100:]
    for name, seq in c_terminals.items():
        alignments = subprocessing.run_muscle_single(name, seq, cterm_file)
        query_seq = alignments[name]
        ref_seq = alignments["EryAII_ref"]
        c_terminal_residues[name] = utils.extract_by_reference_positions(query_seq, ref_seq, [55, 64])
    return c_terminal_residues


def find_possible_orders(genes: List[CDSFeature], start_gene: CDSFeature,
                         end_gene: CDSFeature) -> List[List[CDSFeature]]:
    """ Finds all possible arrangements of the given genes. If not None, the
        start gene will always be the first in each order. Similarly, the end
        gene will always be last.

        Arguments:
            genes: a list of all genes, may include start and end
            start_gene: None or the gene with which to start every arrangement
            end_gene: None or the gene with which to end every arrangement

        Returns:
            a list of lists, each sublist being a unique ordering of the
            provided genes
    """
    assert genes
    genes_to_order = []
    for gene in genes:
        if gene == start_gene or gene == end_gene:
            pass
        else:
            genes_to_order.append(gene)
    possible_orders = []
    start = []  # type: List[CDSFeature]
    if start_gene:
        start = [start_gene]
    end = []  # type: List[CDSFeature]
    if end_gene:
        end = [end_gene]
    for order in list(itertools.permutations(genes_to_order, len(genes_to_order))):
        possible_orders.append(start + list(order) + end)
    return possible_orders


def rank_biosynthetic_orders(n_terminal_residues, c_terminal_residues,
                             possible_orders: List[List[CDSFeature]]) -> List[CDSFeature]:
    """ Scores each possible order according to terminal pairs of adjacent genes.

        Arguments:
            n_terminal_residues: a dictionary mapping genes to their pair of N terminal residues
            c_terminal_residues: a dictionary mapping genes to their pair of C terminal residues
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
    best_order = None
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


def perform_docking_domain_analysis(genes: List[CDSFeature]) -> List[CDSFeature]:
    """ Estimates gene ordering based on docking domains of features

        Arguments:
            genes: a list of genes to order

        Returns:
            a list of genes in estimated order
    """
    start_gene, end_gene = find_first_and_last_genes(genes)
    data_dir = path.get_full_path(__file__, "data", "terminals")

    n_terminal_residues = extract_nterminus(data_dir, genes, start_gene)
    c_terminal_residues = extract_cterminus(data_dir, genes, end_gene)
    possible_orders = find_possible_orders(genes, start_gene, end_gene)

    geneorder = rank_biosynthetic_orders(n_terminal_residues, c_terminal_residues, possible_orders)
    return geneorder


def find_colinear_order(genes: List[CDSFeature]) -> List[CDSFeature]:
    """ Estimates gene ordering based on colinearity

        Arguments:
            genes: a list of genes to order

        Returns:
            a list of genes in estimated order
    """
    direction = 0
    for gene in genes:
        direction += gene.strand
    geneorder = list(genes)
    # Reverse if first gene encodes a multidomain protein with a TE/TD domain
    if direction < 0:
        geneorder.reverse()
    gene_domains = geneorder[0].nrps_pks.domain_names
    if "Thioesterase" in gene_domains or "TD" in gene_domains:
        if len(gene_domains) > 1:
            geneorder.reverse()
    return geneorder
