# License: GNU Affero General Public License v3 or later
# A copy of GNU AGPL v3 should have been included in this software package in LICENSE.txt.

""" A collection of functions for reinterpreting/converting the combined
    predictions into a final, potentially modified, result.
"""

from collections import defaultdict
import logging
import string
from typing import Dict, List, Tuple

from antismash.common.secmet.features import CDSFeature

from .data_structures import Prediction
from .pks_names import get_short_form

ALLOWABLE_PREDICTION_CHARACTERS = set(string.ascii_letters + string.digits).union(set('-,()'))

AVAILABLE_SMILES_PARTS = {'GLY', 'ALA', 'VAL', 'LEU', 'ILE', 'MET', 'PRO', 'PHE', 'TRP', 'SER', 'THR', 'ASN', 'GLN',
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


def calculate_individual_consensus(predictions: List[str]) -> str:
    """ Finds the most frequent prediction in the list of predictions that is
        also in the set of available smiles parts.

        Defaults to 'pk' in the case of no valid predictions or a tie.

        Arguments:
            predictions: the list of predictions

        Returns:
            the most frequent prediction
    """
    best = "pk"
    highest_count = -1
    for pred in set(predictions):
        count = predictions.count(pred)
        if pred not in AVAILABLE_SMILES_PARTS:
            continue
        if count > 0 and count > highest_count:
            best = pred
            highest_count = count
        elif count == highest_count:  # a tie means we default back to pk
            best = "pk"
    return best


def generate_nrps_consensus(results: Dict[str, Prediction]) -> str:
    """ Finds the most frequent prediction in the list of predictions that is
        also in the set of available smiles parts.

        Defaults to 'nrp' in the case of no valid predictions or a tie.

        Arguments:
            predictions: the list of predictions

        Returns:
            the most frequent prediction
    """
    assert isinstance(results, dict), type(results)
    hit_counts = defaultdict(int)  # type: Dict[str, int]
    for method, prediction in results.items():
        assert isinstance(method, str), method
        assert isinstance(prediction, Prediction), prediction
        if len(prediction.get_classification()) == 1:
            best = prediction.get_classification()[0]
            if not set(best).issubset(ALLOWABLE_PREDICTION_CHARACTERS):
                raise ValueError("%s generated bad prediction string: %s" % (method, best))
            hit_counts[best] += 1

    consensus = "nrp"
    # only really care about the first two for checking if it was a tie
    best_hits = sorted((count, name) for name, count in hit_counts.items())[:2]
    # if the best hit isn't tie, use that
    if best_hits and len({count for count, _ in best_hits}) == 1:
        consensus = best_hits[0][1]
    return consensus


def calculate_consensus_prediction(cds_features: List[CDSFeature], results: Dict[str, Dict[str, Prediction]]
                                   ) -> Tuple[Dict[str, str], Dict[str, str]]:
    """ Uses all calculations to generate smiles parts to use

        Arguments:
            cds_features: a list of CDSFeature to calculate consensus for
            results: a dictionary mapping AntismashDomain name to
                        a dictionary mapping analysis method to a Prediction

        Returns:
            a tuple of dicts mapping domain label to prediction
                the first for all domains except AT domains in trans-AT clusters
                the second for AT domains in trans-AT clusters
    """
    # Combine substrate specificity predictions into consensus prediction
    cis_at = {}  # type: Dict[str, str]  # feature name -> prediction
    trans_at = {}  # type: Dict[str, str]  # feature name -> prediction

    for cds in cds_features:
        assert cds.region, "Orphaned CDS found: %s" % cds
        for domain in cds.nrps_pks.domains:
            predictions = results[domain.feature_name]
            if 'OTHER' in domain.label:
                continue
            if domain.name == "PKS_AT":
                preds = []  # type: List[str]
                preds.extend(map(get_short_form, predictions["minowa_at"].get_classification()))
                preds.extend(map(get_short_form, predictions["signature"].get_classification()))
                consensus = calculate_individual_consensus(preds)

                if 'transatpks' not in cds.region.products:
                    cis_at[domain.feature_name] = consensus
                else:
                    trans_at[domain.feature_name] = consensus

            if 'transatpks' in cds.region.products and domain.name == "PKS_KS":
                # For chemical display purpose for chemicals from trans-AT PKS gene cluster
                # mal is always assumed for trans-AT
                cis_at[domain.feature_name] = "mal"

            if domain.name in ["AMP-binding", "A-OX"]:
                cis_at[domain.feature_name] = generate_nrps_consensus(predictions)
            elif domain.name == "CAL_domain":
                preds = predictions["minowa_cal"].get_classification()
                pred = get_short_form(preds[0])
                assert isinstance(pred, str)
                if pred in AVAILABLE_SMILES_PARTS:
                    cis_at[domain.feature_name] = pred
                else:
                    logging.debug("missing %s from SMILES parts for domain %s", pred, domain.feature_name)
                    cis_at[domain.feature_name] = "pk"

    return cis_at, trans_at


def find_duplicate_position(domains: List[str], item: str) -> List[int]:
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
            preds: a dict mapping domain label (e.g. nrpspksdomains_SCO123_AT1)
                   to a prediction for that domain
            target: "PKS_KS" or "PKS_AT" for checking AT vs trans-AT
            target_list: a list of positions in the gene's domains where target is found
            lists: a list of lists of positions for KR, DH and ER domains
            mappings: a list of dictionaries mapping a prediction to an altered prediction

        Returns:
            None
    """
    assert len(lists) == len(mappings)
    for idx, target_element in enumerate(target_list):
        key = "nrpspksdomains_{}_{}.{}".format(locus, target, idx + 1)
        for sublist, mapping in zip(lists, mappings):
            for position in sublist:
                if not target_element < position:
                    continue
                if idx + 1 <= len(target_list) - 1 and not position < target_list[idx + 1]:
                    continue
                current = preds[key]
                preds[key] = mapping.get(current, current)


def modify_monomer_predictions(cds_features: List[CDSFeature], predictions: Dict[str, str]) -> None:
    """ Modifies monomer predictions based on domain construction chain. Changes
        the predictions in place.

        Arguments:
            cds_features: a dictionary mapping CDS name to a list of domain names
                          in the order they are found in the CDS
            predictions: a dict mapping domain label (e.g. nrpspksdomains_SCO123_AT1)
                         to a prediction for that domain

        Returns:
            None
    """
    # for modifications, e.g. mal -> ohmal
    # must be the same length and ordering as the lists generation below
    # i.e. for each mapping there must be exactly one domain position list
    mappings = [{"mal": "ohmal", "mmal": "ohmmal", "mxmal": "ohmxmal", "emal": "ohemal"},  # KR
                {"ohmal": "ccmal", "ohmmal": "ccmmal", "ohmxmal": "ccmxmal", "ohemal": "ccemal"},  # DH
                {"ccmal": "redmal", "ccmmal": "redmmal", "ccmxmal": "redmxmal", "ccemal": "redemal"}]  # ER

    for cds in cds_features:
        assert cds.region, "Orphaned CDS found"
        domain_names = cds.nrps_pks.domain_names
        lists = [find_duplicate_position(domain_names, 'PKS_KR'),
                 find_duplicate_position(domain_names, 'PKS_DH'),
                 find_duplicate_position(domain_names, 'PKS_ER')]

        if 'transatpks' not in cds.region.products:
            label = "PKS_AT"
            data = find_duplicate_position(domain_names, 'PKS_AT')
        else:
            label = "PKS_KS"
            data = find_duplicate_position(domain_names, 'PKS_KS')
        # TODO: should the transat predictions be used if relevant?
        update_prediction(cds.get_name(), predictions, label, data, lists, mappings)
