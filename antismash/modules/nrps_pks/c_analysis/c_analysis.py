# License: GNU Affero General Public License v3 or later
# A copy of GNU AGPL v3 should have been included in this software package in LICENSE.txt.


"""
Predict activity and active site motifs for domains in the condensation
domain superfamily (Condensation, Heterocyclization, Epimerization)
"""

import logging

from antismash.common import brawn, path, utils
from antismash.common.secmet import CDSFeature
from antismash.modules.nrps_pks.data_structures import SimplePrediction, Prediction

# Condensation domain alignment from Rausch et al. 2007: 10.1186/1471-2148-7-78
DATA_DIR = path.get_full_path(__file__, "data")
C_DOMAINS_FILENAME = "Rausch_C_203.fa"


def extract_c_domains(cds_features: list[CDSFeature]) -> list[tuple[str, str, str]]:
    """ Returns information for each C superfamily domain found in the record

        Argument:
            cds_features: A list of CDSFeatures

        Returns:
            a list of tuples, one for each C/E domain, containing the feature name,
                    the specific type of domain, and the protein sequence

    """
    results = []
    for cds in cds_features:
        for domain in cds.nrps_pks.domains:
            if domain.name in ["Cglyc", "Condensation_DCL", "Condensation_LCL",
                               "Condensation_Starter", "Condensation_Dual",
                               "Heterocyclization", "Epimerization"]:
                seq = str(cds.translation)[domain.start:domain.end]
                results.append((domain.feature_name, domain.name, seq))
    return results


def is_active(signature: str, c_type: str) -> str:
    """ Predicts whether a C/E/Cy domain is active based on its signature

        Arguments:
            signature: a string of the active site, length 7
            c_type: a string of the specific domain type
    """
    if signature[1] == "H":
        return "active"
    if c_type == "Heterocyclization" and signature[0] == "D":
        return "active"
    return "inactive"



def run_c_analysis(c_domains: list[tuple[str, str, str]]) -> dict[str, dict[str, Prediction]]:
    """ Extract activity and active site signatures from C superfamily domains

        Arguments:
            c_domains: a list of tuples, one for each C/E domain, containing
            the feature name, the specific type of domain, and the protein sequence

        Returns:
            a dictionary mapping
                analysis method name to a dictionary mapping
                    domain name to Prediction
    """
    method_results: dict[str, dict[str, Prediction]] = {}
    activesites: dict[str, Prediction] = {}
    activities: dict[str, Prediction] = {}

    alignment = brawn.get_cached_alignment(C_DOMAINS_FILENAME, DATA_DIR)
    positions = list(range(132, 139))
    for name, c_type, seq in c_domains:
        ref_id = "Bacisubti.NP_389716.1.ppsA_1_glu_starter"
        aligned, ref_aligned = brawn.get_aligned_pair(seq, ref_id, alignment)

        c_signature = utils.extract_by_reference_positions(aligned, ref_aligned, positions)
        if not c_signature:
            activities[name] = SimplePrediction("c_activity", "unknown")
            continue
        activesites[name] = SimplePrediction("c_activesite", c_signature)

        # Check activity
        activity = is_active(c_signature, c_type)
        activities[name] = SimplePrediction("c_activity", activity)

    method_results["c_activity"] = activities
    method_results["c_activesite"] = activesites
    return method_results


def run_nrps_c_predictions(cds_features: list[CDSFeature]) -> dict[str, dict[str, Prediction]]:
    """ Extract C domains and run active site analysis on the given CDS features

        Arguments:
            cds_features: a list of CDSFeatures

        Returns:
            a dictionary mapping
                analysis method name to a dictionary mapping
                    domain name to Prediction
    """
    logging.info("Predicting NRPS C domain activity using active site motifs")
    c_domains = extract_c_domains(cds_features)
    method_results: dict[str, dict[str, Prediction]] = {}
    if c_domains:
        method_results = run_c_analysis(c_domains)
    return method_results
