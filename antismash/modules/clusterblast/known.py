# License: GNU Affero General Public License v3 or later
# A copy of GNU AGPL v3 should have been included in this software package in LICENSE.txt.

import logging
import os
from helperlibs.wrappers.io import TemporaryDirectory

import antismash.common.deprecated as utils
import antismash.common.path as path
from .core import ClusterResult, GeneralResults, parse_all_clusters, \
                  load_clusterblast_database, internal_homology_blast, \
                  create_blast_inputs, run_diamond, \
                  write_raw_clusterblastoutput, score_clusterblast_output

# Tuple is ( binary_name, optional)
_required_binaries = [
    ('blastp', False),
    ('makeblastdb', False),
    ('diamond', False),
]

_required_files = [
    ('knownclusterprots.fasta', False),
    ('knownclusterprots.dmnd', False),
    ('knownclusters.txt', False)
]

def _get_datafile_path(filename):
    data_dir = path.get_full_path(__file__, 'data')
    return os.path.join(data_dir, 'known', filename)

def check_known_prereqs(options):
    "Check if all required applications are around"
    failure_messages = []
    for binary_name, optional in _required_binaries:
        if path.locate_executable(binary_name) is None and not optional:
            failure_messages.append("Failed to locate file: %r" % binary_name)

    for file_name, optional in _required_files:
        if path.locate_file(_get_datafile_path(file_name)) is None and not optional:
            failure_messages.append("Failed to locate file: %r" % file_name)

    return failure_messages

def run_knownclusterblast_on_record(seq_record, options):
    logging.info('Running known cluster search')
    clusters, proteins = load_clusterblast_database(seq_record, searchtype="knownclusterblast")
    if not options.cb_general:
        seq_record.internalhomologygroupsdict = internal_homology_blast(seq_record)
    results = perform_knownclusterblast(options, seq_record, clusters, proteins)
    utils.CODE_SKIP_WARNING()
    #prepare_data(seq_record, options, searchtype="knownclusters")
    utils.CODE_SKIP_WARNING()
    #generate_Storage_for_cb(options, seq_record, searchtype="KnownClusterBlastData")
    return results

def perform_knownclusterblast(options, seq_record, clusters, proteins):
    # Run BLAST on gene cluster proteins of each cluster and parse output
    logging.info("Running DIAMOND knowncluster searches..")
    results = GeneralResults(seq_record.id, search_type="knownclusterblast")

    all_names, all_seqs = [], []
    for cluster in seq_record.get_clusters():
        names, seqs = create_blast_inputs(cluster)
        all_names.extend(names)
        all_seqs.extend(seqs)

    with TemporaryDirectory(change=True) as tempdir:
        utils.writefasta([qcname.replace(" ", "_") for qcname in all_names],
                         all_seqs, "input.fasta")
        run_diamond("input.fasta", _get_datafile_path('knownclusterprots'),
                    tempdir, options)
        with open("input.out", 'r') as fh:
            blastoutput = fh.read()
        write_raw_clusterblastoutput(options.output_dir, blastoutput,
                                     search_type="knownclusterblast")
    minseqcoverage = 40
    minpercidentity = 45
    clusters_by_number, _ = parse_all_clusters(blastoutput, minseqcoverage,
                                              minpercidentity, seq_record)

    allcoregenes = seq_record.get_cds_features()
    for genecluster in seq_record.get_clusters():
        clusternumber = genecluster.get_cluster_number()
        cluster_names_to_queries = clusters_by_number.get(clusternumber, {})
        ranking = score_clusterblast_output(clusters, allcoregenes, cluster_names_to_queries)
        # store results
        cluster_result = ClusterResult(genecluster, ranking)
        results.add_cluster_result(cluster_result, clusters, proteins)

        utils.CODE_SKIP_WARNING()
#        write_clusterblast_output(options, seq_record, knownclusterblastStorage, searchtype="knownclusters")
    utils.CODE_SKIP_WARNING()
#    mibig_protein_homology(blastoutput, seq_record, clusters, options)
#    logging.critical("mibig homology results discarded")
    return results


def mibig_protein_homology(blastoutput, seq_record, clusters, options):
    logging.critical("mibig homology still only writing to file")
    minseqcoverage = 20
    minpercidentity = 20
    _, queries_by_cluster = parse_all_clusters(blastoutput, minseqcoverage,
                                               minpercidentity, seq_record)
    for genecluster in seq_record.get_clusters():
        cluster_number = genecluster.get_cluster_number()
        queries = queries_by_cluster.get(cluster_number, {})

        # Since the BLAST query was only for proteins in the cluster just need to iterate through the keys and generate
        # a file for each of the keys
        outputfolder = os.path.join(options.output_dir, 'knownclusterblast',
                                    "cluster{}".format(cluster_number))
        if not os.path.exists(outputfolder):
            os.mkdir(outputfolder)
        for cluster_protein in queries.values():
            protein_name = cluster_protein.id
            with open(outputfolder + os.sep + protein_name + '_mibig_hits.txt', 'w') as outfile:
                outfile.write('#Protein\tDescription\tMiBIG Cluster\tMiBIG Product'
                              '\tPercent ID\tPercent Coverage\tBLAST Score\t Evalue\n')
                for subject in cluster_protein.subjects.values():
                    gene_id = subject.locus_tag
                    gene_descr = subject.annotation
                    mibig_cluster = subject.genecluster
                    mibig_product = clusters[mibig_cluster][1]
                    percent_id = str(subject.perc_ident)
                    blast_score = str(subject.blastscore)
                    percent_cvg = str(subject.perc_coverage)
                    e_value = str(subject.evalue)
                    outfile.write(gene_id + '\t' + gene_descr + '\t' + mibig_cluster
                                  + '\t' + mibig_product + '\t' + percent_id + '\t' + percent_cvg
                                  + '\t' + blast_score + '\t' + e_value + '\n')
