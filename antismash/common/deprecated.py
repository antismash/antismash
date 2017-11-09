# License: GNU Affero General Public License v3 or later
# A copy of GNU AGPL v3 should have been included in this software package in LICENSE.txt.

"""
This file will be removed as soon as all modules from antiSMASH 4 have been
converted
"""

import logging

# temporary code skip logging # TODO
import inspect
import linecache

import Bio
from Bio.Seq import Seq
# pylint: disable=unused-import
from Bio.SeqFeature import SeqFeature, FeatureLocation # for others importing
from Bio.SeqRecord import SeqRecord
# pylint: enable=unused-import

from antismash.config import get_config
from antismash.common.all_orfs import scan_orfs, sort_orfs
from antismash.common.secmet import CDSFeature, Feature

# pylint: disable=unused-import
from .utils import generate_unique_id, RobustProteinAnalysis
# pylint: enable=unused-import


def CODE_SKIP_WARNING():
    prev = inspect.currentframe().f_back
    logging.critical("skipping code: " + prev.f_code.co_name +"():" \
            + linecache.getline(prev.f_code.co_filename, prev.f_lineno + 1).replace('%', '%%').rstrip())
# end temp


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
    out_file = open(filename, "w")
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
    raise RuntimeError("get_aa_sequence(cds) called, did you mean cds.get_aa_sequence()?")

def get_feature_dict_protein_id(_record):
    raise RuntimeError("get_feature_dict_protein_id(record) called, did you mean record.get_cds_accession_mapping()?")

def get_feature_dict(_seq_record):
    raise RuntimeError("get_feature_dict(record) called, did you mean record.get_cds_name_mapping()?")

def sortdictkeysbyvaluesrev(_data):
    raise RuntimeError("sortdictkeysbyvaluesrev(data) called, did you mean [i[0] for i in sorted(data.items(), key=lambda x: (x[1], x[0]), reverse=True)]?")

def features_overlap(_a, _b):
    raise RuntimeError("features_overlap(a, b) called, did you mean a.overlaps_with(b)?")
