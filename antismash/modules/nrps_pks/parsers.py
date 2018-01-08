# License: GNU Affero General Public License v3 or later
# A copy of GNU AGPL v3 should have been included in this software package in LICENSE.txt.

import logging

from typing import Dict, List, Set, Tuple

LONG_TO_SHORT = {'Malonyl-CoA': 'mal', 'Methylmalonyl-CoA': 'mmal', 'Methoxymalonyl-CoA': 'mxmal',
                 'Ethylmalonyl-CoA': 'emal', 'Isobutyryl-CoA': 'isobut', '2-Methylbutyryl-CoA': '2metbut',
                 'trans-1,2-CPDA': 'trans-1,2-CPDA', 'Acetyl-CoA': 'Acetyl-CoA', 'Benzoyl-CoA': 'benz',
                 'Propionyl-CoA': 'prop', '3-Methylbutyryl-CoA': '3metbut',
                 'CE-Malonyl-CoA': 'cemal', '2-Rhyd-Malonyl-CoA': '2Rhydmal', 'CHC-CoA': 'CHC-CoA',
                 'inactive': 'inactive'}
SHORT_TO_LONG = {val: key for key, val in LONG_TO_SHORT.items()}


def calculate_individual_consensus(predictions: List[str], available_smiles_parts: Set[str]) -> str:
    """ Finds the most frequent prediction in the list of predictions that is
        also in the set of available smiles parts.

        In the case of a tie, the first prediction in the list with the highest
        score is used.

        Arguments:
            predictions: the list of predictions
            available_smiles_parts: the set of predictions which are 'valid'

        Returns:
            the most frequent prediction
    """
    best = "pk"
    highest_count = -1
    for pred in predictions:
        count = predictions.count(pred)
        if count > 1 and count > highest_count and pred in available_smiles_parts:
            best = pred
            highest_count = count
    return best


def calculate_consensus_prediction(genes, results) -> Tuple[Dict[str, str], Dict[str, str]]:
    """ Uses all calculations to generate smiles parts to use

        Arguments:
            genes: a list of CDSFeature to calculate consensus for
            results: a dictionary mapping PKS analysis method to the prediction
                     for that method

        Returns:
            a tuple of dicts mapping domain label to prediction
                the first for AT domains
                the second for KS domains
    """
    # Combine substrate specificity predictions into consensus prediction
    non_trans_at = {}
    trans_at = {}
    available_smiles_parts = {'GLY', 'ALA', 'VAL', 'LEU', 'ILE', 'MET', 'PRO', 'PHE', 'TRP', 'SER', 'THR', 'ASN', 'GLN',
                              'TYR', 'CYS', 'LYS', 'ARG',
                              'HIS', 'ASP', 'GLU', 'MPRO', 'ORN', 'PGLY', 'DAB', 'BALA', 'AEO', 'DHA', 'PIP', 'BMT',
                              'gly', 'ala', 'val', 'leu', 'ile', 'met', 'pro', 'phe', 'trp', 'ser',
                              'thr', 'asn', 'gln', 'tyr', 'cys', 'lys', 'arg', 'his', 'asp', 'glu', 'aaa', 'mpro',
                              'dhb', '2hiva', 'orn', 'pgly', 'dab', 'bala', 'aeo', '4mha', 'pico', 'phg',
                              'dha', 'scy', 'pip', 'bmt', 'adds', 'aad', 'abu', 'hiv', 'dhpg', 'bht', '3-me-glu',
                              '4pPro', 'ala-b', 'ala-d', 'dht', 'Sal', 'tcl', 'lys-b', 'hpg', 'hyv-d',
                              'iva', 'vol', 'mal', 'mmal', 'ohmal', 'redmal', 'mxmal', 'emal', 'nrp', 'pk', 'Gly',
                              'Ala', 'Val', 'Leu', 'Ile', 'Met', 'Pro', 'Phe', 'Trp', 'Ser', 'Thr', 'Asn', 'Gln', 'Tyr',
                              'Cys', 'Lys', 'Arg', 'His', 'Asp', 'Glu', 'Mpro', '23Dhb', '34Dhb', '2Hiva', 'Orn',
                              'Pgly', 'Dab', 'Bala', 'Aeo', '4Mha', 'Pico', 'Aaa', 'Dha', 'Scy', 'Pip',
                              'Bmt', 'Adds', 'DHpg', 'DHB', 'nrp', 'pk'}

    for feature in genes:
        locus = feature.get_name()

        for domain in feature.nrps_pks.domains:
            if 'OTHER' in domain.label:
                continue
            domain_name = locus + domain.label
            if 'transatpks' not in feature.cluster.products:
                if domain.name == "PKS_AT":
                    preds = []
                    if results["minowa_at"].get(domain_name):
                        pred = results["minowa_at"][domain_name][0][0]
                        preds.append(LONG_TO_SHORT.get(pred))
                    if results["signature"].get(domain_name):
                        preds.append(results["signature"][domain_name][0].name.rsplit("_", 1)[-1])
                    consensus = calculate_individual_consensus(preds, available_smiles_parts)
                    non_trans_at[domain_name] = consensus
            else:
                if domain.name == "PKS_AT":
                    preds = []
                    if results["minowa_at"].get(domain_name):
                        pred = results["minowa_at"][domain_name][0][0]
                        preds.append(LONG_TO_SHORT.get(pred))
                    if results["signature"].get(domain_name):
                        preds.append(results["signature"][domain_name][0].name.rsplit("_", 1)[-1])
                    consensus = calculate_individual_consensus(preds, available_smiles_parts)
                    trans_at[domain_name] = consensus
                # For chemical display purpose for chemicals from trans-AT PKS gene cluster
                # mal is always assumed for trans-AT
                elif domain.name == "PKS_KS":
                    non_trans_at[domain_name] = "mal"
            if domain.name in ["AMP-binding", "A-OX"]:
                non_trans_at[domain_name] = "nrp"
            elif domain.name == "CAL_domain":
                pred = results["minowa_cal"][domain_name][0][0]
                pred = LONG_TO_SHORT.get(pred, pred)
                if pred in available_smiles_parts:
                    non_trans_at[domain_name] = pred
                else:
                    logging.critical("missing %s from SMILES parts for domain %s", pred, domain_name)
                    non_trans_at[domain_name] = "pk"
    return non_trans_at, trans_at


def find_duplicate_position(domains, item) -> List[int]:
    """ Finds all indices of elements of domains that are equal to item

        Arguments:
            domains: the sequence of domains to look over
            item: the item to find the positions of

        Returns:
            a list containing the indices of `item` within `domains`
    """
    start_at = -1
    locs = []
    while start_at < len(domains):
        try:
            loc = domains.index(item, start_at + 1)
        except ValueError:
            break
        else:
            locs.append(loc)
            start_at = loc
    return locs


def update_prediction(locus: str, preds: Dict[str, str], target: str,
                      target_list: List[int], lists: List[List[int]],
                      mappings: List[Dict[str, str]]) -> None:
    """ Updates predictions of a gene's domains. Modifies in place.

        Arguments:
            locus: the name of the gene
            preds: a dict mapping domain label (e.g. SCO123_AT1) to prediction
                   for that domain
            target: "_KS" or "_AT" for checking AT vs trans-AT
            target_list: a list of positions in the gene's domains where target is found
            lists: a list of lists of positions for KR, DH and ER domains
            mappings: a list of dictionaries mapping a prediction to an altered prediction

        Returns:
            None
    """
    assert len(lists) == len(mappings)
    for idx, target_element in enumerate(target_list):
        key = locus + target + str(idx + 1)
        for sublist, mapping in zip(lists, mappings):
            for position in sublist:
                if not target_element < position:
                    continue
                if idx + 1 <= len(target_list) - 1 and not position < target_list[idx + 1]:
                    continue
                current = preds[key]
                preds[key] = mapping.get(current, current)


def modify_monomer_predictions(genes, predictions) -> None:
    """ Modifies monomer predictions based on domain construction chain. Changes
        the predictions in place.

        Arguments:
            gene_domains: a dictionary mapping gene name to a list of domain names
                          in the order they are found in the gene
            predictions:

        Returns:
            None

    """

    # for modifications, e.g. mal -> ohmal
    # must be the same length and ordering as the lists generation below
    # i.e. for each mapping there must be exactly one domain position list
    mappings = [{"mal": "ohmal", "mmal": "ohmmal", "mxmal": "ohmxmal", "emal": "ohemal"},  # KR
                {"ohmal": "ccmal", "ohmmal": "ccmmal", "ohmxmal": "ccmxmal", "ohemal": "ccemal"},  # DH
                {"ccmal": "redmal", "ccmmal": "redmmal", "ccmxmal": "redmxmal", "ccemal": "redemal"}]  # ER

    for gene in genes:
        domain_names = gene.nrps_pks.domain_names
        lists = [find_duplicate_position(domain_names, 'PKS_KR'),
                 find_duplicate_position(domain_names, 'PKS_DH'),
                 find_duplicate_position(domain_names, 'PKS_ER')]

        if 'transatpks' not in gene.cluster.products:
            label = "_AT"
            data = find_duplicate_position(domain_names, 'PKS_AT')
        else:
            label = "_KS"
            data = find_duplicate_position(domain_names, 'PKS_KS')
        # TODO: should the transat predictions be used if relevant?
        update_prediction(gene.get_name(), predictions, label, data, lists, mappings)
