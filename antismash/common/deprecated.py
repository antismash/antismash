# License: GNU Affero General Public License v3 or later
# A copy of GNU AGPL v3 should have been included in this software package in LICENSE.txt.

"""
This file will be removed as soon as all modules from antiSMASH 4 have been
converted
"""

# since this code has a limited lifetime, in this file kill any pylint messages
# pylint: disable=missing-docstring,invalid-name

import logging

# temporary code skip logging
import inspect
import linecache

# pylint: disable=unused-import
from Bio.SeqFeature import SeqFeature, FeatureLocation # for others importing
from Bio.SeqRecord import SeqRecord
# pylint: enable=unused-import

from antismash.config import get_config


def CODE_SKIP_WARNING():
    prev = inspect.currentframe().f_back
    logging.critical("skipping code: " + prev.f_code.co_name +"():" \
            + linecache.getline(prev.f_code.co_filename, prev.f_lineno + 1).replace('%', '%%').rstrip())
# end temp


def get_pksnrps_cds_features(seq_record) -> list:
    features = seq_record.get_cds_features_within_clusters()
    pksnrpscoregenes = []
    for feature in features:
        if feature.nrps_pks.domains:
            pksnrpscoregenes.append(feature)
    return pksnrpscoregenes

def get_nrpspks_domain_dict(seq_record) -> dict:
    domaindict = {}
    for feature in seq_record.get_cds_features():
        if feature.nrps_pks.domains:
            domaindict[feature.get_name()] = list(feature.nrps_pks.domains)
    return domaindict


def get_version() -> str:
    return get_config().version


def hmmlengths(hmmfile) -> dict:
    lengths = {}
    with open(hmmfile, "r") as handle:
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


def get_cluster_features_of_type(record, product):
    "Return all cluster features within a record that have a product type"
    return [cluster for cluster in record.get_clusters() if product in cluster.products]


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
    raise RuntimeError("get_aa_sequence(cds) called, did you mean cds.translation?")

def get_feature_dict_protein_id(_record):
    raise RuntimeError("get_feature_dict_protein_id(record) called, did you mean record.get_cds_accession_mapping()?")

def get_feature_dict(_seq_record):
    raise RuntimeError("get_feature_dict(record) called, did you mean record.get_cds_name_mapping()?")

def sortdictkeysbyvaluesrev(_data):
    raise RuntimeError("sortdictkeysbyvaluesrev(data) called, did you mean [i[0] for i in sorted(data.items(), key=lambda x: (x[1], x[0]), reverse=True)]?")

def features_overlap(_a, _b):
    raise RuntimeError("features_overlap(a, b) called, did you mean a.overlaps_with(b)?")
