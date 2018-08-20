# License: GNU Affero General Public License v3 or later
# A copy of GNU AGPL v3 should have been included in this software package in LICENSE.txt.

# since this code has a limited lifetime, in this file kill any pylint messages
# pylint: disable=missing-docstring,invalid-name

"""
This file will be removed as soon as all modules from antiSMASH 4 have been
converted
"""

from typing import no_type_check

# pylint: disable=unused-import
from Bio.SeqFeature import SeqFeature, FeatureLocation  # for others importing
from Bio.SeqRecord import SeqRecord
# pylint: enable=unused-import

# DEAD FUNCTIONS
# these only exist so that the mapping to new functions is easier to do


@no_type_check
def get_pksnrps_cds_features(_seq_record) -> list:
    raise NotImplementedError("get_pksnrps_cds_features(record) called, use record.get_nrps_pks_cds_features()")


@no_type_check
def get_nrpspks_domain_dict(_seq_record) -> dict:
    raise RuntimeError("get_nrpspks_domain_dict() should not be called")


@no_type_check
def get_version() -> str:
    raise NotImplementedError("get_version() called, use AntismashConfig.version")


@no_type_check
def hmmlengths(_hmmfile) -> dict:
    raise NotImplementedError("hmmlengths() called, use antismash.common.utils.get_hmm_lengths(hmmfile)")


@no_type_check
def get_all_features_of_type(_seq_record, _types):
    raise RuntimeError("get_all_features_of_type(record, types) called, did you mean record.get_*()")


@no_type_check
def get_cds_features_within_clusters(_seq_record):
    raise RuntimeError("get_withincluster_cds_features(record) called, use record.get_cds_features_within_regions()")


@no_type_check
def get_withincluster_cds_features(_seq_record):
    raise RuntimeError("get_withincluster_cds_features(record) called, use record.get_cds_features_within_regions()")


@no_type_check
def get_gene_id(_feature):
    raise RuntimeError("using get_gene_id(feature), did you mean feature.get_name() or feature.unique_id")


@no_type_check
def get_aa_translation(_seq_record, _feature):
    raise RuntimeError("get_aa_translation(record, feature) called, use record.get_aa_translation(feature)")


@no_type_check
def get_cluster_type(_cluster):
    raise RuntimeError("get_cluster_type(cluster) called,"
                       " did you mean cluster.get_product_string() or cluster.products?")


@no_type_check
def get_cluster_by_nr(_seq_record, _queryclusternr):
    raise RuntimeError("get_cluster_type(seq_record, cluster_num) called,"
                       " did you mean seq_record.get_cluster(cluster_num)")


@no_type_check
def sort_features(_seq_record):
    raise RuntimeError("utils.sort_features(seq_record) called, did you mean sorted(seq_record.get_all_features())?")


@no_type_check
def get_cluster_cds_features(_cluster, _seq_record):
    raise RuntimeError("utils.get_cluster_cds_features(cluster) called, did you mean cluster.cds_children?")


@no_type_check
def get_aa_sequence(_feature, **_kwargs):
    raise RuntimeError("get_aa_sequence(cds) called, did you mean cds.translation?")


@no_type_check
def get_feature_dict_protein_id(_record):
    raise RuntimeError("get_feature_dict_protein_id(record) called, did you mean record.get_cds_accession_mapping()?")


@no_type_check
def get_feature_dict(_seq_record):
    raise RuntimeError("get_feature_dict(record) called, did you mean record.get_cds_name_mapping()?")


@no_type_check
def sortdictkeysbyvaluesrev(_data):
    raise RuntimeError("sortdictkeysbyvaluesrev(data) called, did you mean "
                       "[i[0] for i in sorted(data.items(), key=lambda x: (x[1], x[0]), reverse=True)]?")


@no_type_check
def features_overlap(_a, _b):
    raise RuntimeError("features_overlap(a, b) called, did you mean a.overlaps_with(b)?")
