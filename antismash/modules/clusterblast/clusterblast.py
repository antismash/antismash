# License: GNU Affero General Public License v3 or later
# A copy of GNU AGPL v3 should have been included in this software package in LICENSE.txt.

import logging
import os

from helperlibs.wrappers.io import TemporaryDirectory

from .core import write_fastas_with_all_genes, run_diamond, write_raw_clusterblastoutput, \
                  parse_all_clusters, score_clusterblast_output
from .results import ClusterResult, GeneralResults

def perform_clusterblast(options, record, db_clusters, db_proteins) -> GeneralResults:
    """ Run BLAST on gene cluster proteins for each cluster, parse output and
        return result rankings for each cluster

        Arguments:
            options: antismash Config
            record: the Record to analyse
            db_clusters: a dict mapping reference cluster name to ReferenceCluster
            db_proteins: a dict mapping reference protein to Protein

        Returns:
            a GeneralResults instance with results for each cluster in the record
    """
    clusters = record.get_clusters()
    with TemporaryDirectory(change=True) as tempdir:
        write_fastas_with_all_genes(clusters, "input.fasta")
        run_diamond("input.fasta",
                    os.path.join(options.database_dir, 'clusterblast', 'geneclusterprots'),
                    tempdir, options)

        with open("input.out", 'r') as handle:
            blastoutput = handle.read()

        write_raw_clusterblastoutput(options.output_dir, blastoutput)

        clusters_by_number, _ = parse_all_clusters(blastoutput, record,
                                                   min_seq_coverage=10,
                                                   min_perc_identity=30)
        results = GeneralResults(record.id)

        for cluster in clusters:
            cluster_number = cluster.get_cluster_number()
            cluster_names_to_queries = clusters_by_number.get(cluster_number, {})
            allcoregenes = [cds.get_accession() for cds in record.get_cds_features()]
            ranking = score_clusterblast_output(db_clusters, allcoregenes, cluster_names_to_queries)

            # store the results
            result = ClusterResult(cluster, ranking, db_proteins)
            results.add_cluster_result(result, db_clusters, db_proteins)

        results.write_to_file(record, options)
    return results
