# License: GNU Affero General Public License v3 or later
# A copy of GNU AGPL v3 should have been included in this software package in LICENSE.txt.

""" Calculates a likely order of NRPS/PKS domains """

import itertools
import logging
from typing import Dict, List, Optional, Tuple

from antismash.common import brawn, path, utils
from antismash.common.secmet import CDSFeature, Module, Record
from antismash.detection.nrps_pks_domains.modular_domain import ModularDomain

from .html_output import will_handle
from .results import CandidateClusterPrediction, modify_substrate
from .smiles_generator import gen_smiles_from_pksnrps

C_TERMINAL_PATH = path.get_full_path(__file__, "data", "terminals", "cterm.fasta")
N_TERMINAL_PATH = path.get_full_path(__file__, "data", "terminals", "nterm.fasta")


def find_split_module_chains(cds_features: List[CDSFeature], record: Record) -> Dict[CDSFeature, CDSFeature]:
    """ Builds a mapping of module head to module tail for each cross-CDS module
        in the given features

        Arguments:
            cds_features: the CDS features to find cross-CDS modules in

        Returns:
            a dictionary mapping each cross-CDS module head CDS to the tail CDS
    """
    chains: Dict[CDSFeature, CDSFeature] = {}  # track ordering of CDS features forming a complete module
    for cds in cds_features:
        for module in cds.modules:
            if len(module.parent_cds_names) > 1:
                parents = (record.get_cds_by_name(name) for name in module.parent_cds_names)
                head, tail = sorted(parents, reverse=cds.location.strand == -1)
                chains[head] = tail
    return chains


def get_follower_genes(cdses: list[CDSFeature]) -> Dict[str, str]:
    """ Returns a mapping of gene name to immediately subsequent gene name, if
        on the same strand. If no valid follower is present for a gene, the gene
        won't be present in the mapping.

        Arguments:
            cdses: a list of CDS features

        Returns:
            a dictionary mapping gene name to following gene name
    """
    neighbours = {}
    for i, cds in enumerate(cdses):
        if cds.location.strand == -1:
            next_index = i - 1
            if next_index < 0:
                continue
        else:
            next_index = i + 1
            if next_index >= len(cdses):
                continue
        next_cds = cdses[next_index]
        if next_cds.location.strand == cds.location.strand:
            neighbours[cds.get_name()] = next_cds.get_name()
    return neighbours


def analyse_biosynthetic_order(nrps_pks_features: List[CDSFeature],
                               consensus_predictions: Dict[str, str],
                               record: Record) -> List[CandidateClusterPrediction]:
    """ For each NRPS or PKS candidate cluster, determines if that candidate cluster is
        docking or not then determines the monomer ordering

        Arguments:
            nrps_pks_features: all NRPS/PKS features within the record
            consensus_predictions: a dictionary mapping each NRPS/PKS domain name to its prediction
            record: the Record being analysed

        Returns:
            a dictionary mapping candidate cluster number to
                a tuple of
                    prediction string
                    and whether docking domain analysis was used for the prediction
    """
    compound_predictions: List[CandidateClusterPrediction] = []
    # Find NRPS/PKS gene candidate_clusters
    candidate_clusters = [cluster for cluster in record.get_candidate_clusters()
                          if will_handle(cluster.products, cluster.product_categories)]
    if not candidate_clusters:
        return []
    # Predict biosynthetic gene order in candidate clusters using starter domains,
    # thioesterase domains, gene order and docking domains
    for candidate_cluster in candidate_clusters:
        candidate_cluster_number = candidate_cluster.get_candidate_cluster_number()
        cds_in_candidate_cluster = [gene for gene in nrps_pks_features if gene.overlaps_with(candidate_cluster)]
        if not cds_in_candidate_cluster:
            continue
        pks_features, nrps_count, hybrid_count = find_candidate_cluster_modular_enzymes(cds_in_candidate_cluster)
        pks_chains = {}
        if not hybrid_count:
            pks_chains = find_split_module_chains(pks_features, record)
        # use docking domain analysis, if possible, to identify order
        # since this will grow as n!, an upper limit is required
        if 2 <= len(pks_features) < 11 + len(pks_chains) and not nrps_count and not hybrid_count:
            logging.debug("CandidateCluster %d monomer ordering method: domain docking analysis",
                          candidate_cluster_number)
            neighbours = get_follower_genes(cds_in_candidate_cluster)
            geneorder = perform_docking_domain_analysis(pks_features, pks_chains, neighbours)
            docking = True
        else:
            logging.debug("CandidateCluster %d monomer ordering method: colinear", candidate_cluster_number)
            with_complete = filter(lambda cds: any(module.is_complete() for module in cds.modules),
                                   cds_in_candidate_cluster)
            geneorder = find_colinear_order(list(with_complete))
            docking = False

        polymer, smiles = generate_substrates_order(geneorder, consensus_predictions)
        gene_names_in_order = [cds.get_name() for cds in geneorder]
        prediction = CandidateClusterPrediction(candidate_cluster_number, polymer,
                                                docking, smiles, gene_names_in_order)
        compound_predictions.append(prediction)
    return compound_predictions


def find_candidate_cluster_modular_enzymes(cds_features: List[CDSFeature]) -> Tuple[List[CDSFeature], int, int]:
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
    def is_pks_module(module: Module) -> bool:
        domain_names = set(domain.domain for domain in module.domains)
        return bool(domain_names.intersection({"PKS_KS", "PKS_AT"}))

    def is_nrps_module(module: Module) -> bool:
        domain_names = set(domain.domain for domain in module.domains)
        return bool(domain_names.intersection({"AMP-binding", "A-OX", "Condensation"}))

    pks_features = []
    nrps_count = 0
    hybrid_count = 0
    for cds in cds_features:
        has_pks = any(is_pks_module(module) for module in cds.modules if module.is_complete())
        has_nrps = any(is_nrps_module(module) for module in cds.modules if module.is_complete())
        if has_nrps and has_pks:
            hybrid_count += 1
        elif has_pks:
            pks_features.append(cds)
        elif has_nrps:
            nrps_count += 1
    return pks_features, nrps_count, hybrid_count


def generate_substrates_order(geneorder: List[CDSFeature], consensus_predictions: Dict[str, str]
                              ) -> Tuple[str, str]:
    """ Generate substrates order and SMILES from predicted gene order and consensus
        predictions. E.g. (ala-dpg) + (pk).

        Arguments:
            geneorder: a list of CDSFeatures
            consensus_predictions: a dictionary mapping domain name to prediction

        Returns:
            a tuple of the polymer and the smiles for the polymer, both as strings
    """
    components = []
    monomers_by_cds = []

    for gene in geneorder:
        monomers = []
        for module in gene.modules:
            if not module.is_complete():
                continue
            # for cross-CDS modules, report the monomer only for the head CDS
            if len(module.parent_cds_names) > 1 and module.parent_cds_names[0] != gene.get_name():
                continue
            substrate = ""
            for domain in module.domains:
                consensus = consensus_predictions.get(domain.get_name())
                if consensus:
                    substrate = consensus
                    break
            monomer = modify_substrate(module, substrate)
            if not monomer:
                continue
            monomers.append(monomer)
            domain_names = []
            for domain in module.domains:
                assert isinstance(domain, ModularDomain)
                if domain.domain == "MT" and domain.domain_subtype:
                    domain_names.append(domain.domain_subtype)
                    continue
                domain_names.append(domain.domain or "")
            components.append((substrate, monomer, domain_names))

        if monomers:
            monomers_by_cds.append(f"({' - '.join(monomers)})")

    polymer = " + ".join(monomers_by_cds)
    smiles = gen_smiles_from_pksnrps(components)

    return polymer, smiles


def find_first_and_last_cds(cds_features: List[CDSFeature]) -> Tuple[Optional[CDSFeature], Optional[CDSFeature]]:
    """ Find first and last CDSFeature based on the CDS having starter or
        finalisation modules.

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
        if cds.modules and cds.modules[-1].is_final_module():
            module = cds.modules[-1]
            # if this CDS contains the head of a cross-CDS module, it can't be the end
            if len(module.parent_cds_names) > 1 and cds.get_name() == module.parent_cds_names[0]:
                continue
            # two possible ends, this really ought to be two products
            if end_cds:
                end_cds = None
                break
            end_cds = cds

    # find the start
    for cds in cds_features:
        if cds == end_cds:
            continue
        if cds.modules and cds.modules[0].is_starter_module():
            if start_cds:
                # two possible starts, don't attempt fallbacks
                return None, end_cds
            start_cds = cds

    # if no starter module found, try looking for KS-AT-ACP
    if not start_cds:
        for cds in cds_features:
            if cds == end_cds:
                continue
            domain_names = cds.nrps_pks.domain_names
            if domain_names[:3] in [["PKS_KS", "PKS_AT", "ACP"], ["PKS_KS", "PKS_AT", "PKS_PP"]]:
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
    for cds in cds_features:
        if cds is not start_cds:
            seq = str(cds.translation)
            n_terminals[cds.get_name()] = seq[:50]
    reference_name = "EryAIII_5_6_ref"
    alignment = brawn.get_cached_alignment(N_TERMINAL_PATH, data_dir)
    for name, seq in n_terminals.items():
        query_seq, ref_seq = brawn.get_aligned_pair(seq, reference_name, alignment)
        residue = utils.extract_by_reference_positions(query_seq, ref_seq, [2, 15])
        assert residue
        n_terminal_residues[name] = residue
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
    c_terminals: Dict[str, str] = {}
    for cds in cds_features:
        if cds is not end_cds:
            seq = str(cds.translation)
            c_terminals[cds.get_name()] = seq[-100:]
    reference_name = "EryAII_ref"
    alignment = brawn.get_cached_alignment(C_TERMINAL_PATH, data_dir)
    for name, seq in c_terminals.items():
        query_seq, ref_seq = brawn.get_aligned_pair(seq, reference_name, alignment)
        residue = utils.extract_by_reference_positions(query_seq, ref_seq, [55, 64])
        assert residue
        c_terminal_residues[name] = residue
    return c_terminal_residues


def find_possible_orders(cds_features: List[CDSFeature], start_cds: Optional[CDSFeature],
                         end_cds: Optional[CDSFeature], chains: Dict[CDSFeature, CDSFeature]
                         ) -> List[List[CDSFeature]]:
    """ Finds all possible arrangements of the given cds_features. If not None, the
        start gene will always be the first in each order. Similarly, the end
        gene will always be last.

        Arguments:
            cds_features: a list of all CDSFeatures, may include start_cds and end_cds
            start_cds: None or the CDS with which to start every arrangement
            end_cds: None or the CDS with which to end every arrangement
            chains: a dictionary mapping cross-CDS module head CDS to tail CDS

        Returns:
            a list of lists, each sublist being a unique ordering of the
            provided CDSFeatures
    """
    assert len(cds_features) - len(chains) < 11, "input too large, function is O(n!)"
    assert start_cds is None or isinstance(start_cds, CDSFeature)
    assert end_cds is None or isinstance(end_cds, CDSFeature)
    if start_cds or end_cds:
        assert start_cds != end_cds, "Using same gene for start and end of ordering"

    tails = {v: k for k, v in chains.items()}

    # if the end CDS is the tail end of a split module,
    # use the beginning of that chain as the 'end CDS'
    if end_cds:
        while end_cds in tails:
            end_cds = tails[end_cds]

    cds_to_order = []
    for cds in cds_features:
        # the permutations shouldn't include chained CDSes or the start or end
        if cds in (start_cds, end_cds) or cds in tails:
            continue
        cds_to_order.append(cds)

    possible_orders = []
    start: List[CDSFeature] = []
    if start_cds:
        start = [start_cds]
    end: List[CDSFeature] = []
    if end_cds:
        # if the end CDS is the tail end of a split module,
        # use the beginning of that chain as the 'end CDS'
        while end_cds in tails:
            end_cds = tails[end_cds]
        end = [end_cds]
    for order in itertools.permutations(cds_to_order, len(cds_to_order)):
        # expand out any subchains formed by cross-CDS modules
        full = []
        for cds in start + list(order) + end:
            full.append(cds)
            next_cds = chains.get(cds)
            while next_cds:
                full.append(next_cds)
                next_cds = chains.get(next_cds)
        possible_orders.append(full)
    # ensure the list of possible orders is itself ordered for reliability
    return sorted(possible_orders, key=lambda x: [g.location.start for g in x])


def rank_biosynthetic_orders(n_terminal_residues: Dict[str, str],
                             c_terminal_residues: Dict[str, str],
                             possible_orders: List[List[CDSFeature]],
                             following_genes: Dict[str, str]) -> List[CDSFeature]:
    """ Scores each possible order according to terminal pairs of adjacent cds_features.

        Arguments:
            n_terminal_residues: a dictionary mapping CDSFeature to their pair of N terminal residues
            c_terminal_residues: a dictionary mapping CDSFeature to their pair of C terminal residues
            possible_orders: a list of gene orderings to evaluate
            following_genes: a mapping of gene name to subsequent immediate neighbouring
                             gene name if on the same strand

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
            if following_genes.get(gene.get_name()) == next_gene.get_name():
                score += 1
            # additional CDS features may be brought in if they formed part of a
            # cross-CDS module and contained no other complete modules
            if gene.get_name() not in c_terminal_residues or next_gene.get_name() not in n_terminal_residues:
                continue
            res1a, res2a = tuple(c_terminal_residues[gene.get_name()])
            res1b, res2b = tuple(n_terminal_residues[next_gene.get_name()])
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


def perform_docking_domain_analysis(cds_features: List[CDSFeature], chains: Dict[CDSFeature, CDSFeature],
                                    following_genes: dict[str, str]) -> List[CDSFeature]:
    """ Estimates gene ordering based on docking domains of features

        Arguments:
            cds_features: a list of CDSFeatures to order
            chains: a dictionary mapping cross-CDS module head CDS to tail CDS
            following_genes: a mapping of gene name to subsequent immediate neighbouring
                             gene name if on the same strand

        Returns:
            a list of CDSFeatures in estimated order
    """
    start_cds, end_cds = find_first_and_last_cds(cds_features)
    data_dir = path.get_full_path(__file__, "data")

    n_terminal_residues = extract_nterminus(data_dir, cds_features, start_cds)
    c_terminal_residues = extract_cterminus(data_dir, cds_features, end_cds)
    possible_orders = find_possible_orders(cds_features, start_cds, end_cds, chains)

    geneorder = rank_biosynthetic_orders(n_terminal_residues, c_terminal_residues,
                                         possible_orders, following_genes)
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
    if not geneorder:
        return geneorder
    # Reverse if first gene encodes a multidomain protein with a TE/TD domain
    if direction < 0:
        geneorder.reverse()
    gene_domains = geneorder[0].nrps_pks.domain_names
    if "Thioesterase" in gene_domains or "TD" in gene_domains:
        if len(gene_domains) > 1:
            geneorder.reverse()
    return geneorder
