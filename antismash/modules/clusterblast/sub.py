# License: GNU Affero General Public License v3 or later
# A copy of GNU AGPL v3 should have been included in this software package in LICENSE.txt.

import functools
import logging
from multiprocessing import Process, Pool
import os

from helperlibs.wrappers.io import TemporaryDirectory

from antismash.common.deprecated import CODE_SKIP_WARNING, writefasta
import antismash.common.path as path
from .core import write_raw_clusterblastoutput, \
                  create_blast_inputs, blastparse, score_clusterblast_output, \
                  write_clusterblast_output, read_clusterblast_output, runblast, \
                  load_clusterblast_database, internal_homology_blast, \
                  ClusterResult, GeneralResults

def _get_datafile_path(filename):
    data_dir = path.get_full_path(__file__, 'data')
    return os.path.join(data_dir, 'sub', filename)

def check_sub_prereqs(options):
    "Check if all required applications and datafiles are around"
    # Tuple is ( binary_name, optional)
    _required_binaries = [
        ('blastp', False),
        ('makeblastdb', False),
    ]

    _required_files = [
        ('subclusterprots.fasta', False),
        ('subclusterprots.fasta.phr', False),
        ('subclusterprots.fasta.pin', False),
        ('subclusterprots.fasta.psq', False),
        ('subclusters.txt', False)
    ]
    failure_messages = []
    for binary_name, optional in _required_binaries:
        if path.locate_executable(binary_name) is None and not optional:
            failure_messages.append("Failed to locate file: %r" % binary_name)

    for file_name, optional in _required_files:
        if path.locate_file(_get_datafile_path(file_name)) is None and not optional:
            failure_messages.append("Failed to locate file: %r" % file_name)

    return failure_messages

def run_subclusterblast_on_record(seq_record, options):
    logging.info('Running subcluster search')
    clusters, proteins = load_clusterblast_database(seq_record, "subclusterblast")
    if not options.cb_general:
        seq_record.internalhomologygroupsdict = internal_homology_blast(seq_record)
    results = perform_subclusterblast(options, seq_record, clusters, proteins)
    CODE_SKIP_WARNING()
    #prepare_data(seq_record, options, searchtype="subclusterblast")
    CODE_SKIP_WARNING()
    #generate_Storage_for_cb(options, seq_record, searchtype="SubClusterBlastData")
    return results


def write_clusterblast_inputfiles(options, queryclusternames, queryclusterseqs):
    size = int(len(queryclusternames) / options.cpus)
    for i in range(options.cpus):
        if i == 0:
            setnames = queryclusternames[:size]
            setseqs = queryclusterseqs[:size]
        elif i == (options.cpus - 1):
            setnames = queryclusternames[i * size:]
            setseqs = queryclusterseqs[i * size:]
        else:
            setnames = queryclusternames[i * size:(i + 1) * size]
            setseqs = queryclusterseqs[i * size:(i + 1) * size]
        writefasta(setnames, setseqs, "input" + str(i) + ".fasta")

def _run_blast_helper(database, index):
    """A wrapper of runblast() just to reverse args for functools.partial use"""
    # note: can't be nested because it must be picklable for Pool's use
    return runblast("input%d.fasta" % index, database)

def run_clusterblast_processes(options):
    database = _get_datafile_path('subclusterprots.fasta')
    # set the first arg to always be database
    partial = functools.partial(_run_blast_helper, database)
    # then run over each file in parallel
    with Pool(options.cpus) as pool:
        pool.map(partial, range(options.cpus))

def perform_subclusterblast(options, seq_record, clusters, proteins):
    #Run BLAST on gene cluster proteins of each cluster and parse output
    logging.info("Running NCBI BLAST+ subcluster searches..")
    results = GeneralResults(seq_record.id, search_type="subclusterblast")
    with TemporaryDirectory(change=True):
        allcoregenes = [cds.get_accession() for cds in seq_record.get_cds_features()]
        for cluster in seq_record.get_clusters():
            clusternumber = cluster.get_cluster_number()
            logging.info("   Gene cluster " + str(clusternumber))
            # prepare and run diamond
            queryclusternames, queryclusterseqs = create_blast_inputs(cluster)
            write_clusterblast_inputfiles(options, queryclusternames, queryclusterseqs)
            run_clusterblast_processes(options)
            blastoutput = read_clusterblast_output(options)
            write_raw_clusterblastoutput(options.output_dir, blastoutput, search_type="subclusterblast")
            logging.info("   Blast search finished. Parsing results...")

            # parse and score diamond results
            minseqcoverage = 40
            minpercidentity = 45
            _, cluster_names_to_queries = blastparse(blastoutput, minseqcoverage, minpercidentity, seq_record)
            ranking = score_clusterblast_output(clusters, allcoregenes, cluster_names_to_queries)
            logging.critical("cluster %s has %d subclusterblast results", cluster.location, len(ranking))
            # store results
            cluster_result = ClusterResult(cluster, ranking)
            results.add_cluster_result(cluster_result, clusters, proteins)
            CODE_SKIP_WARNING()
            #write_clusterblast_output(options, seq_record, subclusterblastStorage, searchtype="subclusterblast")
    return results
