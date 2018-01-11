# License: GNU Affero General Public License v3 or later
# A copy of GNU AGPL v3 should have been included in this software package in LICENSE.txt.

"""
    Runs minowa and kr_streochem over PKS domains
"""

from antismash.common.secmet import CDSFeature

from .results import PKSResults
from .substrates_pks import count_pks_genes, run_minowa_predictor_pks_at, \
                            run_minowa_predictor_pks_cal, run_kr_stereochemistry_predictions, \
                            extract_at_domains


def run_pks_substr_spec_predictions(genes: CDSFeature) -> PKSResults:
    """ Runs the PKS analyses on the given genes and returns the various predictions
        as a PKSResults object
    """
    at_domains = extract_at_domains(genes)
    counted = count_pks_genes(genes)
    results = PKSResults()
    if at_domains:
        signature_results, minowa_at_results = run_minowa_predictor_pks_at(at_domains)
        results.method_results["signature"] = signature_results
        results.method_results["minowa_at"] = minowa_at_results
    if counted:
        minowa_cal_results = run_minowa_predictor_pks_cal(genes)
        kr_activity, kr_stereo = run_kr_stereochemistry_predictions(genes)
        results.method_results["minowa_cal"] = minowa_cal_results
        results.method_results["kr_activity"] = kr_activity
        results.method_results["kr_stereochem"] = kr_stereo
    return results
