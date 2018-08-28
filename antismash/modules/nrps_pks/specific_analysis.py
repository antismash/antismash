#  License: GNU Affero General Public License v3 or later
#  A copy of GNU AGPL v3 should have been included in this software package in LICENSE.txt.

'''
In-depth analysis and annotation of NRPS/PKS regions.
'''

import logging
from typing import List

from antismash.common.secmet import Record, CDSFeature, AntismashDomain, Region
from antismash.config import ConfigType

from .orderfinder import analyse_biosynthetic_order
from .parsers import calculate_consensus_prediction, modify_monomer_predictions
from .results import NRPS_PKS_Results
from .smiles_generator import gen_smiles_from_pksnrps
from .substrates import run_pks_substr_spec_predictions

from .nrps_predictor import run_nrpspredictor


def get_a_domains_from_cds_features(record: Record, cds_features: List[CDSFeature]) -> List[AntismashDomain]:
    """ Fetches all AMP-binding AntismashDomains which are contained within the given
        CDS features.

        Arguments:
            record: the Record containing both AntismashDomains and CDSFeatures
            cds_features: the specific CDSFeatures from which to get the A-domains

        Returns:
            a list of AntismashDomains, one for each A domain found
    """
    a_domains = []
    for cds in cds_features:
        for domain in cds.nrps_pks.domains:
            if domain.name in ["AMP-binding", "A-OX"]:
                as_domain = record.get_domain_by_name(domain.feature_name)
                assert isinstance(as_domain, AntismashDomain), type(as_domain)
                a_domains.append(as_domain)
    return a_domains


def specific_analysis(record: Record, results: NRPS_PKS_Results, options: ConfigType) -> NRPS_PKS_Results:
    """ Runs the various NRPS/PKS analyses on a record and returns their results """
    nrps_pks_genes = record.get_nrps_pks_cds_features()

    if not nrps_pks_genes:
        logging.debug("No NRPS or PKS genes found, skipping analysis")
        return results

    a_domains = get_a_domains_from_cds_features(record, nrps_pks_genes)
    if a_domains:
        logging.info("Predicting A domain substrate specificities with NRPSPredictor2")
        results.add_method_results("NRPSPredictor2", run_nrpspredictor(a_domains, options))

    pks_results = run_pks_substr_spec_predictions(nrps_pks_genes)
    for method, method_results in pks_results.items():
        results.add_method_results(method, method_results)
    results.consensus, results.consensus_transat = calculate_consensus_prediction(nrps_pks_genes, results.domain_predictions)

    modify_monomer_predictions(nrps_pks_genes, results.consensus)

    supercluster_predictions = analyse_biosynthetic_order(nrps_pks_genes, results.consensus, record)
    for prediction in supercluster_predictions:
        prediction.smiles = gen_smiles_from_pksnrps(prediction.polymer)
        supercluster = record.get_supercluster(prediction.supercluster_number)
        region = supercluster.parent
        assert isinstance(region, Region), type(region)
        results.region_predictions[region.get_region_number()].append(prediction)
    return results
