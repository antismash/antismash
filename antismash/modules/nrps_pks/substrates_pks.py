# License: GNU Affero General Public License v3 or later
# A copy of GNU AGPL v3 should have been included in this software package in LICENSE.txt.

import logging

from helperlibs.wrappers.io import TemporaryDirectory

from .minowa import minowa_cal, minowa_at
from .kr_analysis import kr_analysis
from .at_analysis import at_analysis


def count_pks_genes(genes):
    """ returns the combined count of PKS_AT, PKS_KR, and CAL_domain """
    pkscount = 0
    for gene in genes:
        for domain in gene.nrps_pks.domains:
            if domain.name in ["PKS_AT", "CAL_domain", "PKS_KR"]:
                pkscount += 1
    return pkscount


def extract_pks_genes(genes):
    """ returns a list of names and sequences for each PKS_AT domain found in the record """
    results = {}
    for gene in genes:
        locus = gene.get_name()
        domains = gene.nrps_pks.domains
        count = 0
        for domain in domains:
            if domain.name == "PKS_AT":
                count += 1
                seq = str(gene.get_aa_sequence())[domain.start:domain.end]
                name = locus + "_AT" + str(count)
                results[name] = seq
    return results


def run_minowa_predictor_pks_at(pks_domains):
    """ analyses AT domains with Minowa and signature based detection """
    # Predict PKS AT domain specificities with Minowa et al. method and PKS code (NP searcher / ClustScan / own?)
    # Run PKS signature analysis
    logging.info("Predicting PKS AT domain substrate specificities by Yadav et al. PKS signature sequences")
    signature_results = at_analysis.run_at_domain_analysis(pks_domains)

    # Minowa method: run Minowa_AT
    logging.info("Predicting PKS AT domain substrate specificities by Minowa et al. method")
    with TemporaryDirectory(change=True):
        minowa_results = minowa_at.run_minowa_at(pks_domains)
    return signature_results, minowa_results


def run_minowa_predictor_pks_cal(genes):
    """ Predict PKS CAL domain specificities with Minowa et al. method. """
    cal_domains = {}
    logging.info("Predicting CAL domain substrate specificities by Minowa et al. method")
    for gene in genes:
        locus = gene.get_name()
        count = 0
        for domain in gene.nrps_pks.domains:
            if domain.name == "CAL_domain":
                count += 1
                seq = str(gene.get_aa_sequence())[domain.start:domain.end]
                name = locus + "_CAL" + str(count)
                cal_domains[name] = seq
    with TemporaryDirectory(change=True):
        minowa_results = minowa_cal.run_minowa_cal(cal_domains)
    return minowa_results


def run_kr_stereochemistry_predictions(genes):
    """ Predict PKS KR domain stereochemistry using pattern as published in ClustScan
    """
    queries = {}
    logging.info("Predicting PKS KR activity and stereochemistry using KR "
                 "fingerprints from Starcevic et al.")
    for gene in genes:
        locus = gene.get_name()
        count = 0
        for domain in gene.nrps_pks.domains:
            if domain.name == "PKS_KR":
                count += 1
                seq = str(gene.get_aa_sequence())[domain.start:domain.end]
                name = locus + "_KR" + str(count)
                queries[name] = seq
    activity, stereo = kr_analysis.run_kr_analysis(queries)
    return activity, stereo
