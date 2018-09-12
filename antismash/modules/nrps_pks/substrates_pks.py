# License: GNU Affero General Public License v3 or later
# A copy of GNU AGPL v3 should have been included in this software package in LICENSE.txt.

""" The section of nrps_pks responsible for analysing PKS components """

import logging
from typing import Dict, List, Tuple

from helperlibs.wrappers.io import TemporaryDirectory

from antismash.common.secmet import CDSFeature

from .data_structures import Prediction
from .minowa import minowa_cal, minowa_at
from .kr_analysis import kr_analysis
from .at_analysis import at_analysis


def count_pks_genes(genes: List[CDSFeature]) -> int:
    """ returns the combined count of PKS_AT, PKS_KR, and CAL_domain """
    pkscount = 0
    for gene in genes:
        for domain in gene.nrps_pks.domains:
            if domain.name in ["PKS_AT", "CAL_domain", "PKS_KR"]:
                pkscount += 1
    return pkscount


def extract_at_domains(cds_features: List[CDSFeature]) -> Dict[str, str]:
    """ Returns a dictionary mapping domain name to domain sequence
        for each PKS_AT domain found in the record
    """
    results = {}
    for cds in cds_features:
        for domain in cds.nrps_pks.domains:
            if domain.name == "PKS_AT":
                seq = str(cds.translation)[domain.start:domain.end]
                results[domain.feature_name] = seq
    return results


def run_minowa_predictor_pks_at(at_domains: Dict[str, str]
                                ) -> Tuple[Dict[str, Prediction],
                                           Dict[str, Prediction]]:
    """ analyses AT domains with Minowa and signature based detection """
    # Run PKS signature analysis
    logging.info("Predicting PKS AT domain substrate specificities by Yadav et al. PKS signature sequences")
    signature_results = at_analysis.run_at_domain_analysis(at_domains)

    # Minowa method
    logging.info("Predicting PKS AT domain substrate specificities by Minowa et al. method")
    with TemporaryDirectory(change=True):
        minowa_results = minowa_at.run_minowa_at(at_domains)
    return signature_results, minowa_results


def run_minowa_predictor_pks_cal(cds_features: List[CDSFeature]) -> Dict[str, Prediction]:
    """ Predict PKS CAL domain specificities with Minowa et al. method. """
    cal_domains = {}
    logging.info("Predicting CAL domain substrate specificities by Minowa et al. method")
    for cds in cds_features:
        for domain in cds.nrps_pks.domains:
            if domain.name == "CAL_domain":
                seq = str(cds.translation)[domain.start:domain.end]
                cal_domains[domain.feature_name] = seq
    with TemporaryDirectory(change=True):
        results = minowa_cal.run_minowa_cal(cal_domains)
    return results


def run_kr_stereochemistry_predictions(cds_features: List[CDSFeature]) -> Tuple[Dict[str, Prediction],
                                                                                Dict[str, Prediction]]:
    """ Predict PKS KR domain stereochemistry using pattern as published in ClustScan
    """
    queries = {}
    logging.info("Predicting PKS KR activity and stereochemistry using KR "
                 "fingerprints from Starcevic et al.")
    for cds in cds_features:
        for domain in cds.nrps_pks.domains:
            if domain.name == "PKS_KR":
                seq = str(cds.translation)[domain.start:domain.end]
                queries[domain.feature_name] = seq
    activity, stereo = kr_analysis.run_kr_analysis(queries)
    return activity, stereo
