# License: GNU Affero General Public License v3 or later
# A copy of GNU AGPL v3 should have been included in this software package in LICENSE.txt.

""" The general clusterblast, comparing clusters to other antismash-predicted
    clusters.
"""

import os
from typing import Dict

from antismash.config import ConfigType
from antismash.common.secmet import Record

from .core import (
    get_core_gene_ids,
    parse_all_clusters,
    run_diamond_on_all_regions,
    score_clusterblast_output,
)
from .data_structures import ReferenceCluster, Protein
from .results import RegionResult, GeneralResults


def perform_clusterblast(options: ConfigType, record: Record,
                         db_clusters: Dict[str, ReferenceCluster],
                         db_proteins: Dict[str, Protein]) -> GeneralResults:
    """ Run BLAST on gene cluster proteins for each cluster, parse output and
        return result rankings for each cluster

        Arguments:
            options: antismash Config
            record: the Record to analyse
            db_clusters: a dict mapping reference cluster name to ReferenceCluster
            db_proteins: a dict mapping reference protein name to Protein

        Returns:
            a GeneralResults instance with results for each cluster in the record
    """
    regions = record.get_regions()
    database = os.path.join(options.database_dir, 'clusterblast', 'proteins')
    blastoutput = run_diamond_on_all_regions(regions, database)

    clusters_by_number, _ = parse_all_clusters(blastoutput, record,
                                               min_seq_coverage=10,
                                               min_perc_identity=30)
    results = GeneralResults(record.id)

    core_gene_accessions = get_core_gene_ids(record)

    for region in regions:
        region_number = region.get_region_number()
        cluster_names_to_queries = clusters_by_number.get(region_number, {})
        ranking = score_clusterblast_output(db_clusters, core_gene_accessions,
                                            cluster_names_to_queries)

        # store the results
        result = RegionResult(region, ranking, db_proteins, "general")
        results.add_region_result(result, db_clusters, db_proteins)

    results.write_to_file(record, options)
    return results
