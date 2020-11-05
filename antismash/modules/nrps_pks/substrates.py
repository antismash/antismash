# License: GNU Affero General Public License v3 or later
# A copy of GNU AGPL v3 should have been included in this software package in LICENSE.txt.

"""
    Contains a runner for all PKS-based prediction methods
"""

from typing import Dict, List

from antismash.common.secmet import CDSFeature

from .data_structures import Prediction
from .substrates_pks import (
    count_pks_genes,
    extract_at_domains,
    run_minowa_predictor_pks_at,
    run_minowa_predictor_pks_cal,
    run_kr_stereochemistry_predictions,
)


def run_pks_substr_spec_predictions(cds_features: List[CDSFeature]) -> Dict[str, Dict[str, Prediction]]:
    """ Runs all PKS analyses on the given CDS features

        Arguments:
            cds_features: a list of CDSFeature to run predictions on

        Returns:
            a dictionary mapping
                analysis method name to a dictionary mapping
                    domain name to Prediction
    """
    at_domains = extract_at_domains(cds_features)
    method_results: Dict[str, Dict[str, Prediction]] = {}
    if at_domains:
        signature_results, minowa_at_results = run_minowa_predictor_pks_at(at_domains)
        method_results["signature"] = signature_results
        method_results["minowa_at"] = minowa_at_results
    if count_pks_genes(cds_features):
        minowa_cal_results = run_minowa_predictor_pks_cal(cds_features)
        kr_activity, kr_stereo = run_kr_stereochemistry_predictions(cds_features)
        method_results["minowa_cal"] = minowa_cal_results
        method_results["kr_activity"] = kr_activity
        method_results["kr_stereochem"] = kr_stereo
    return method_results
