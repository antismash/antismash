# License: GNU Affero General Public License v3 or later
# A copy of GNU AGPL v3 should have been included in this software package in LICENSE.txt.

import logging
import os

import antismash.common.deprecated as utils
from helperlibs.wrappers.io import TemporaryDirectory

from .core import create_blast_inputs, run_diamond, write_raw_clusterblastoutput, \
                  parse_all_clusters, score_clusterblast_output
from .results import ClusterResult, GeneralResults

def write_fasta_with_all_genes(clusters, filename):
    all_names, all_seqs = [], []
    for cluster in clusters:
        names, seqs = create_blast_inputs(cluster)
        all_names.extend(names)
        all_seqs.extend(seqs)
    utils.writefasta(all_names, all_seqs, filename)

def perform_clusterblast(options, seq_record, db_clusters, db_proteins):
    #Run BLAST on gene cluster proteins of each cluster and parse output
    geneclusters = seq_record.get_clusters()
    with TemporaryDirectory(change=True) as tempdir:
        write_fasta_with_all_genes(geneclusters, "input.fasta")
        run_diamond("input.fasta",
                    os.path.join(options.database_dir, 'clusterblast', 'geneclusterprots'),
                    tempdir, options)

        with open("input.out", 'r') as handle:
            blastoutput = handle.read()

        write_raw_clusterblastoutput(options.output_dir, blastoutput)

        clusters_by_number, _ = parse_all_clusters(blastoutput, seq_record,
                                                   minseqcoverage=10,
                                                   minpercidentity=30)
        results = GeneralResults(seq_record.id)

        for genecluster in geneclusters:
            cluster_number = genecluster.get_cluster_number()
            cluster_names_to_queries = clusters_by_number.get(cluster_number, {})
            allcoregenes = [cds.get_accession() for cds in seq_record.get_cds_features()]
            ranking = score_clusterblast_output(db_clusters, allcoregenes, cluster_names_to_queries)

            # store the results
            result = ClusterResult(genecluster, ranking, db_proteins)
            results.add_cluster_result(result, db_clusters, db_proteins)

        results.write_to_file(seq_record, options)
    return results
