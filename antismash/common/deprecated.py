# License: GNU Affero General Public License v3 or later
# A copy of GNU AGPL v3 should have been included in this software package in LICENSE.txt.

"""
This file will be removed as soon as all modules from antiSMASH 4 have been
converted
"""

import logging
import os
import re
import sys

import Bio

from Bio.Seq import Seq
from Bio.SeqFeature import SeqFeature, FeatureLocation # for others importing
from Bio.SeqRecord import SeqRecord
from helperlibs.bio import seqio

from antismash.common import gff_parser
from antismash.common.all_orfs import scan_orfs, sort_orfs
from antismash.common.secmet import Record, CDSFeature, Feature

from .utils import generate_unique_id, RobustProteinAnalysis

# temporary code skip logging # TODO
import inspect
import linecache

def CODE_SKIP_WARNING():
    prev = inspect.currentframe().f_back
    logging.critical("skipping code:" + prev.f_code.co_name +"():" \
            + linecache.getline(prev.f_code.co_filename, prev.f_lineno + 1).replace('%', '%%'))
# end temp


def get_feature_dict(seq_record) -> dict:
    """Get a dictionary mapping features to their IDs"""
    features = seq_record.get_cds_features()
    feature_by_id = {}
    for feature in features:
        gene_id = feature.get_name()
        feature_by_id[gene_id] = feature
    return feature_by_id


def get_multifasta(seq_record) -> str:
    """Extract multi-protein FASTA from all CDS features in sequence record"""
    features = seq_record.get_cds_features()
    all_fastas = []
    for feature in features:
        gene_id = feature.get_name()
        fasta_seq = feature.translation
        if "-" in str(fasta_seq):
            fasta_seq = Seq(str(fasta_seq).replace("-", ""), Bio.Alphabet.generic_protein)

        # Never write empty fasta entries
        if not fasta_seq:
            logging.error("No translation for CDS %s", gene_id)
            raise ValueError("No translation for CDS %s" % gene_id)

        all_fastas.append(">%s\n%s" % (gene_id, fasta_seq))
    full_fasta = "\n".join(all_fastas)
    return full_fasta

def writefasta(names, seqs, filename) -> None:
    "Write sequence to a file"
    e = 0
    f = len(names) - 1
    out_file = open(filename,"w")
    while e <= f:
        out_file.write(">")
        out_file.write(names[e])
        out_file.write("\n")
        out_file.write(seqs[e])
        out_file.write("\n")
        e += 1
    out_file.close()

def strip_record(seq_record) -> None:
    """ Discard antismash specific features and feature qualifiers """
    seq_record.clear_clusters()
    seq_record.clear_cluster_borders()
    seq_record.clear_cds_motifs()
    seq_record.clear_antismash_domains()

    # clean up antiSMASH annotations in CDS features
    for feature in seq_record.get_cds_features():
        feature.sec_met = None


def get_smcog_annotations(seq_record):
    logging.critical("get_smcog_annotations(): should use secmet for smCOG note")
    smcogdict = {}
    smcogdescriptions = {}
    for feature in seq_record.get_cds_features():
        for note in feature.notes:
            if "smCOG: " in note:
                smcogid = note.partition("smCOG: ")[2].partition(":")[0]
                smcog_descr = note.partition("smCOG: ")[2].partition(":")[2].partition("(Score:")[0]
                smcogdict[feature.get_name()] = smcogid
                smcogdescriptions[smcogid] = smcog_descr
    return smcogdict, smcogdescriptions

def get_pksnrps_cds_features(seq_record) -> list:
    logging.critical("skipping get_pksnrps_cds_features(), missing PKSNRPS module")
    return []
#    features = get_cds_features(seq_record)
#    pksnrpscoregenes = []
#    for feature in features:
#        if 'sec_met' in feature.qualifiers:
#            for annotation in feature.qualifiers['sec_met']:
#                if annotation.startswith('NRPS/PKS Domain:'):
#                    pksnrpscoregenes.append(feature)
#                    break
#    return pksnrpscoregenes

def get_nrpspks_domain_dict(seq_record) -> dict:
    logging.critical("skipping get_pksnrps_domain_dict(), missing PKSNRPS module")
    return {}
#    domaindict = {}
#    features = get_cds_features(seq_record)
#    for feature in features:
#        domainlist = []
#        if 'sec_met' in feature.qualifiers:
#            domains = [qualifier for qualifier in feature.qualifiers['sec_met'] if "NRPS/PKS Domain: " in qualifier]
#            for domain in domains:
#                hit_id =  domain.partition("NRPS/PKS Domain: ")[2].partition(" (")[0]
#                domstart = domain.partition("NRPS/PKS Domain: ")[2].partition(" (")[2].partition("-")[0]
#                domend = domain.partition("NRPS/PKS Domain: ")[2].partition("). ")[0].rpartition("-")[2]
#                evalue = domain.partition("E-value: ")[2].partition(". Score:")[0]
#                bitscore = domain.partition("Score: ")[2].partition(";")[0]
#                domainlist.append([hit_id, int(domstart), int(domend), evalue, float(bitscore)])
#            if len(domainlist) > 0:
#                domaindict[feature.get_name()] = domainlist
#    return domaindict

def get_nrpspks_substr_spec_preds(seq_record):
    logging.critical("skipping get_nrpspks_substr_spec_preds(), missing PKSNRPS module")
    class Dummy:
        pass
    substr_spec_preds = Dummy()
    substr_spec_preds.consensuspreds = {}
    substr_spec_preds.nrps_svm_preds = {}
    substr_spec_preds.nrps_code_preds = {}
    substr_spec_preds.minowa_nrps_preds = {}
    substr_spec_preds.pks_code_preds = {}
    substr_spec_preds.minowa_pks_preds = {}
    substr_spec_preds.minowa_cal_preds = {}
    substr_spec_preds.kr_activity_preds = {}
    substr_spec_preds.kr_stereo_preds = {}
    return substr_spec_preds
#    features = get_cds_features(seq_record)
#    for feature in features:
#        nrat, nra, nrcal, nrkr = 0, 0, 0, 0
#        if 'sec_met' in feature.qualifiers:
#            domains = [qualifier for qualifier in feature.qualifiers['sec_met'] if "NRPS/PKS Domain: " in qualifier]
#            for domain in domains:
#                if "AMP-binding" in domain or "A-OX" in domain:
#                    nra += 1
#                    domainname = feature.get_name() + "_A" + str(nra)
#                    predictionstext = domain.partition("Substrate specificity predictions:")[2]
#                    nrps_svm_pred = predictionstext.partition(" (NRPSPredictor2 SVM)")[0]
#                    nrps_code_pred = predictionstext.partition(" (NRPSPredictor2 SVM), ")[2].partition(" (Stachelhaus code)")[0]
#                    minowa_nrps_pred = predictionstext.partition("(Stachelhaus code), ")[2].partition(" (Minowa)")[0]
#                    consensuspred = predictionstext.partition("(Minowa), ")[2].partition(" (consensus)")[0]
#                    substr_spec_preds.nrps_svm_preds[domainname] = nrps_svm_pred
#                    substr_spec_preds.nrps_code_preds[domainname] = nrps_code_pred
#                    substr_spec_preds.minowa_nrps_preds[domainname] = minowa_nrps_pred
#                    substr_spec_preds.consensuspreds[domainname] = consensuspred
#                elif "PKS_AT" in domain:
#                    nrat += 1
#                    domainname = feature.get_name() + "_AT" + str(nrat)
#                    predictionstext = domain.partition("Substrate specificity predictions:")[2]
#                    pks_code_pred = predictionstext.partition(" (PKS signature)")[0]
#                    minowa_pks_pred = predictionstext.partition("(PKS signature), ")[2].partition(" (Minowa)")[0]
#                    consensuspred = predictionstext.partition("(Minowa), ")[2].partition(" (consensus)")[0]
#                    substr_spec_preds.pks_code_preds[domainname] = pks_code_pred
#                    substr_spec_preds.minowa_pks_preds[domainname] = minowa_pks_pred
#                    substr_spec_preds.consensuspreds[domainname] = consensuspred
#                elif "CAL_domain" in domain:
#                    nrcal += 1
#                    domainname = feature.get_name() + "_CAL" + str(nrcal)
#                    predictionstext = domain.partition("Substrate specificity predictions:")[2]
#                    minowa_cal_pred = predictionstext.partition(" (Minowa)")[0]
#                    consensuspred = predictionstext.partition("(Minowa), ")[2].partition(" (consensus)")[0]
#                    substr_spec_preds.minowa_cal_preds[domainname] = minowa_cal_pred
#                    substr_spec_preds.consensuspreds[domainname] = consensuspred
#                elif "PKS_KR" in domain:
#                    nrkr += 1
#                    domainname = feature.get_name() + "_KR" + str(nrkr)
#                    activityprediction = domain.partition("Predicted KR activity: ")[2].partition(";")[0]
#                    stereoprediction = domain.partition("Predicted KR stereochemistry: ")[2].partition(";")[0]
#                    substr_spec_preds.kr_activity_preds[domainname] = activityprediction
#                    substr_spec_preds.kr_stereo_preds[domainname] = stereoprediction
#    return substr_spec_preds

def get_structure_pred(cluster) -> str:
    "Return all a structure prediction for a cluster feature"
    for note in cluster.notes:
        if "Monomers prediction: " in note:
            return note.partition("Monomers prediction: ")[2]
    if cluster.get_product_string() == 'ectoine':
        return 'ectoine'
    return "N/A"

def get_version() -> str:
    logging.critical("dummy get_version() being called")
    return "antismash-5.alpha"

def find_all_orfs(seq_record, cluster) -> list: # the old lassopeptides.find_all_orfs
    """Find all ORFs in gene cluster outside annotated CDS features"""
    # Get sequence just for the gene cluster
    fasta_seq = seq_record.seq[cluster.location.start:cluster.location.end]

    # Find orfs throughout the cluster
    forward_matches = scan_orfs(fasta_seq, 1, cluster.location.start)
    reverse_matches = scan_orfs(fasta_seq.complement(), -1, cluster.location.start)
    all_orfs = forward_matches + reverse_matches

    orfnr = 1
    new_features = []

    for orf in sort_orfs(all_orfs):
        # Remove if overlaps with existing CDSs
        skip = False
        for cds in cluster.cds_children:
            if orf.start in cds.location or orf.end in cds.location or cds.location.start in orf or cds.location.end in orf:
                skip = True
                break
        if skip:
            continue
        loc = orf
        dummy_feature = Feature(loc, feature_type="dummy")
        feature = CDSFeature(loc, str(seq_record.get_aa_translation_of_feature(dummy_feature)),
                             locus_tag='cluster_%s_allorf%03d' % (cluster.get_cluster_number(), orfnr))
        new_features.append(feature)
        orfnr += 1

    return new_features

def distance_to_pfam(seq_record, query, hmmer_profiles) -> int: #also from lassopeptides
    """Function to check how many nt a gene is away from a gene with one of a list of given Pfams"""
    nt = 40000 #maximum number of nucleotides distance to search
    #Get all CDS features in seq_record
    cds_features = seq_record.get_cds_features()
    #Get all CDS features within <X nt distances
    close_cds_features = []
    distance = {}
    for cds in cds_features:
        if query.location.start - nt <= cds.location.start <= query.location.end + nt or \
           query.location.start - nt <= cds.location.end <= query.location.end + nt:
            close_cds_features.append(cds)
            distance[cds.get_name()] = min([
                                abs(cds.location.start - query.location.end),
                                abs(cds.location.end - query.location.start),
                                abs(cds.location.start - query.location.start),
                                abs(cds.location.end - query.location.end)])
    #For nearby CDS features, check if they have hits to the pHMM
    closest_distance = -1
    for cds in close_cds_features:
        if cds.sec_met:
            for profile in hmmer_profiles:
                if profile in cds.sec_met.domains:
                    if closest_distance == -1 or distance[cds.get_name()] < closest_distance:
                        closest_distance = distance[cds.get_name()]
    return closest_distance

def get_specific_multifasta(features) -> str:
    """Extract multi-protein FASTA from provided features"""
    all_fastas = []
    for feature in features:
        all_fastas.append(">%s\n%s" % (feature.get_name(), feature.translation))
    return "\n".join(all_fastas)

def hmmlengths(hmmfile) -> dict:
    lengths = {}
    with open(hmmfile,"r") as handle:
        contents = handle.read()
    contents = contents.replace("\r", "\n")
    hmms = contents.split("//")[:-1]
    for hmm in hmms:
        namepart = hmm.split("NAME  ")[1]
        name = namepart.split("\n")[0]
        lengthpart = hmm.split("LENG  ")[1]
        length = lengthpart.split("\n")[0]
        lengths[name] = int(length)
    return lengths

# DEAD FUNCTIONS
# these only exist so that the mapping to new functions is easier to do
def get_all_features_of_type(_seq_record, _types):
    raise RuntimeError("get_all_features_of_type(record, types) called, did you mean record.get_*()")

def get_cds_features_within_clusters(_seq_record):
    raise RuntimeError("get_withincluster_cds_features(record) called, use record.get_cds_features_within_clusters()")

def get_withincluster_cds_features(_seq_record):
    raise RuntimeError("get_withincluster_cds_features(record) called, use record.get_cds_features_within_clusters()")

def get_gene_id(_feature):
    raise RuntimeError("using get_gene_id(feature), did you mean feature.get_name() or feature.unique_id")

def get_aa_translation(_seq_record, _feature):
    raise RuntimeError("get_aa_translation(record, feature) called, use record.get_aa_translation(feature)")

def get_cluster_type(_cluster):
    raise RuntimeError("get_cluster_type(cluster) called, did you mean cluster.get_product_string() or cluster.products?")

def get_cluster_by_nr(_seq_record, _queryclusternr):
    raise RuntimeError("get_cluster_type(seq_record, cluster_num) called, did you mean seq_record.get_cluster(cluster_num)")

def sort_features(_seq_record):
    raise RuntimeError("utils.sort_features(seq_record) called, did you mean sorted(seq_record.get_all_features())?")

def get_cluster_cds_features(_cluster, _seq_record):
    raise RuntimeError("utils.get_cluster_cds_features(cluster) called, did you mean cluster.cds_children?")

def get_aa_sequence(_feature, **_kwargs):
    raise RuntimeError("get_aa_sequence(cds) called, did you mean cds.get_aa_sequence()?")

def get_feature_dict_protein_id(_record):
    raise RuntimeError("get_feature_dict_protein_id(record) called, did you mean record.get_cds_mapping()?")

def sortdictkeysbyvaluesrev(_data):
    raise RuntimeError("sortdictkeysbyvaluesrev(data) called, did you mean [i[0] for i in sorted(data.items(), key=lambda x: (x[1], x[0]), reverse=True)]?")

def features_overlap(_a, _b):
    raise RuntimeError("features_overlap(a, b) called, did you mean a.overlaps_with(b)?")
