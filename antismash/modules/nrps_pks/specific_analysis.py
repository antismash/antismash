#  License: GNU Affero General Public License v3 or later
#  A copy of GNU AGPL v3 should have been included in this software package in LICENSE.txt.

'''
In-depth analysis and annotation of NRPS/PKS gene clusters.
'''

import logging

from Bio.SeqFeature import FeatureLocation

from antismash.common import deprecated
from antismash.common.secmet import AntismashDomain

from .orderfinder import analyse_biosynthetic_order
from .parsers import calculate_consensus_prediction, modify_monomer_predictions, LONG_TO_SHORT
from .structure_drawer import generate_chemical_structure_preds
from .substrates import run_pks_substr_spec_predictions

from .results import PKSResults

DOMAIN_TYPE_MAPPING = {'Condensation_DCL': 'Condensation',
                       'Condensation_LCL': 'Condensation',
                       'Condensation_Dual': 'Condensation',
                       'Condensation_Starter': 'Condensation',
                       'CXglyc': 'Condensation',
                       'Cglyc': 'Condensation',
                       'cMT': 'MT',
                       'oMT': 'MT',
                       'nMT': 'MT',
                       'Polyketide_cyc': 'Polyketide_cyc',
                       'Polyketide_cyc2': 'Polyketide_cyc'}


_UNKNOWN = "(unknown)"


def generate_structure_images(record, results, options):
    """ Generate the structure images based on monomers prediction for all
        cluster features
    """
    compound_predictions = {key: val[0] for key, val in results.cluster_predictions.items()}
    if compound_predictions:
        generate_chemical_structure_preds(compound_predictions, record, options)


def specific_analysis(record, results, options):
    nrps_pks_genes = deprecated.get_pksnrps_cds_features(record)

    if not nrps_pks_genes:
        logging.debug("No NRPS or PKS genes found, skipping analysis")
        return results

    logging.critical("no NRPS prediction methods to use")

    results.pks = run_pks_substr_spec_predictions(nrps_pks_genes)
    results.consensus, results.consensus_transat = calculate_consensus_prediction(nrps_pks_genes, results.pks.method_results)

    modify_monomer_predictions(nrps_pks_genes, results.consensus)

    results.cluster_predictions = analyse_biosynthetic_order(nrps_pks_genes, results.consensus, record)
    generate_structure_images(record, results, options)
    write_data_to_record(nrps_pks_genes, record, results)
    return results


def write_data_to_record(nrps_pks_genes, record, results):
    """ Save substrate specificity predictions in NRPS/PKS domain sec_met info of record

        Workaround to extract positional information for CDS_motifs from the sec_met qualifiers
    """

    logging.critical("the awful nrps specific_analysis.write_data_to_record() is running")

    for feature in nrps_pks_genes:
        at_count = 0
        a_count = 0
        cal_count = 0
        kr_count = 0
        x_count = 0
        nrps_qualifier = feature.nrps_pks
        new_features = []
        gene_id = feature.get_name()
        for domain in nrps_qualifier.domains:
            domain_type = domain.name
            start_aa = domain.start
            end_aa = domain.end
            evalue = domain.evalue
            score = domain.bitscore

            domain.predictions.clear()

            # calculate respective positions based on aa coordinates
            if feature.location.strand == 1:
                start = feature.location.start + (3 * start_aa)
                end = feature.location.start + (3* end_aa)
            else:
                end = feature.location.end - (3 * start_aa)
                start = feature.location.end - (3 * end_aa)
            loc = FeatureLocation(start, end, strand=feature.strand)

            #  set up new CDS_motif feature
            new_feature = AntismashDomain(loc)
            new_feature.domain_subtype = domain_type
            if feature.locus_tag:
                new_feature.locus_tag = feature.locus_tag
            else:
                new_feature.locus_tag = gene_id
            new_feature.detection = "hmmscan"
            new_feature.database = "nrpspksdomains.hmm"
            new_feature.evalue = evalue
            new_feature.score = score
            if feature.transl_table:
                transl_table = feature.transl_table
            else:
                transl_table = 1
            new_feature.translation = str(new_feature.extract(record.seq).translate(table=transl_table))
            if domain_type == "AMP-binding":
                a_count += 1
                domainname = gene_id + "_A" + str(a_count)
                new_feature.label = domainname
                new_feature.domain_id = "nrpspksdomains_" + domainname
                domain.predictions["consensus"] = "nrp"

            elif domain_type == "PKS_AT":
                at_count += 1
                domainname = gene_id + "_AT" + str(at_count)
                new_feature.label = domainname
                new_feature.domain_id = "nrpspksdomains_" + domainname

                # For t1pks, t2pks and t3pks
                if 'transatpks' not in feature.cluster.products:
                    consensus = results.consensus[domainname]
                else: # for transatpks
                    consensus = results.consensus_transat[domainname]
                pks_sig = results.pks.method_results["signature"][domainname]
                if pks_sig:
                    domain.predictions["PKS signature"] = pks_sig[0].name.rsplit("_", 1)[1]
                else:
                    domain.predictions["PKS signature"] = _UNKNOWN
                minowa = results.pks.method_results["minowa_at"][domainname][0][0]
                domain.predictions["Minowa"] = LONG_TO_SHORT.get(minowa, minowa)
                domain.predictions["consensus"] = consensus

            elif domain_type == "CAL_domain":
                cal_count += 1
                domainname = gene_id + "_CAL" + str(cal_count)
                new_feature.label = domainname
                new_feature.domain_id = "nrpspksdomains_" + domainname
                minowa = results.pks.method_results["minowa_cal"][domainname][0][0]
                domain.predictions["Minowa"] = LONG_TO_SHORT.get(minowa, minowa)

            elif domain_type == "PKS_KR":
                kr_count += 1
                domainname = gene_id + "_KR" + str(kr_count)
                new_feature.label = domainname
                new_feature.domain_id = "nrpspksdomains_" + domainname

                domain.predictions["KR activity"] = "active" if results.pks.method_results["kr_activity"][domainname] else "inactive"
                domain.predictions["KR stereochemistry"] = results.pks.method_results["kr_stereochem"].get(domainname, _UNKNOWN)
            else:
                x_count += 1
                new_feature.domain_id = "nrpspksdomains_" + gene_id.partition(".")[0] + "_Xdom"+'{:02d}'.format(x_count)
#                updated_nrps_qualifier.append(domain) # TODO weird, but should it be done?
            for method, pred in domain.predictions.items():
                new_feature.specificity.append("%s: %s" % (method, pred))
            mapping = DOMAIN_TYPE_MAPPING.get(domain_type)
            if mapping:
                new_feature.domain_subtype = domain_type
                new_feature.domain = mapping
            new_features.append(new_feature)

        for new_feature in new_features:
            record.add_feature(new_feature)
