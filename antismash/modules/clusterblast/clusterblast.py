# License: GNU Affero General Public License v3 or later
# A copy of GNU AGPL v3 should have been included in this software package in LICENSE.txt.

import logging
import os

import antismash.common.deprecated as utils
from helperlibs.wrappers.io import TemporaryDirectory

from .core import create_blast_inputs, run_diamond, write_raw_clusterblastoutput, \
                  parse_all_clusters, score_clusterblast_output
from .results import ClusterResult, GeneralResults

def perform_clusterblast(options, seq_record, db_clusters, db_proteins):
    #Run BLAST on gene cluster proteins of each cluster and parse output
    geneclusters = seq_record.get_clusters()
    with TemporaryDirectory(change=True) as tempdir:
        all_names, all_seqs = [], []
        for genecluster in geneclusters:
            names, seqs = create_blast_inputs(genecluster)
            all_names.extend(names)
            all_seqs.extend(seqs)
        logging.info("Running DIAMOND gene cluster search..")
        utils.writefasta(all_names, all_seqs, "input.fasta")
        diamond_result = run_diamond("input.fasta",
                 os.path.join(options.database_dir, 'clusterblast', 'geneclusterprots'),
                 tempdir, options)
        logging.info("   DIAMOND search finished. Parsing results...")

        with open("input.out", 'r') as handle:
            blastoutput = handle.read()

        write_raw_clusterblastoutput(options.output_dir, blastoutput)


        minseqcoverage = 10
        minpercidentity = 30
        clusters_by_number, _ = parse_all_clusters(blastoutput, minseqcoverage,
                                                   minpercidentity, seq_record)
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
