# License: GNU Affero General Public License v3 or later
# A copy of GNU AGPL v3 should have been included in this software package in LICENSE.txt.

""" The knownclusterblast variant of clusterblast, comparing clusters to MIBiG
    clusters.
"""

import logging
import os
from typing import Dict, List

from antismash.common.path import find_latest_database_version
from antismash.common.secmet import Record
from antismash.config import ConfigType, get_config

from .core import (
    check_clusterblast_files,
    get_core_gene_ids,
    load_clusterblast_database,
    parse_all_clusters,
    run_diamond_on_all_regions,
    score_clusterblast_output,
)
from .results import RegionResult, GeneralResults, write_clusterblast_output
from .data_structures import MibigEntry, ReferenceCluster, Protein


def _get_datafile_path(filename: str, config: ConfigType) -> str:
    """ A helper to construct absolute paths to files in the knownclusterblast
        data directory.

        Arguments:
            filename: the name only of the file

        Returns:
            the absolute path of the file
    """
    root = os.path.join(config.database_dir, "knownclusterblast")
    version = find_latest_database_version(root)
    return os.path.join(root, version, filename)


def check_known_prereqs(options: ConfigType) -> List[str]:
    """ Determines if any prerequisite data files or executables are missing

        Arguments:
            options: antismash Config

        Returns:
            a list of error messages, one for each failing prequisite check
    """
    failure_messages = []
    for binary_name in ['blastp', 'makeblastdb', 'diamond']:
        if binary_name not in options.executables:
            failure_messages.append(f"Failed to locate file: {binary_name!r})")

    return failure_messages


def prepare_known_data(logging_only: bool = False) -> List[str]:
    """ Prepare diamond database for KnownClusterBlast.

        Arguments:
            logging_omly: Only log errors, don't try to regenerate databases

        Returns:
            a list of error messages, one for each failing prerequisite check
    """
    failure_messages = []
    config = get_config()
    if "mounted_at_runtime" in config.database_dir:
        return []
    cluster_defs = _get_datafile_path("clusters.txt", config)
    protein_seqs = _get_datafile_path("proteins.fasta", config)
    db_file = _get_datafile_path("proteins.dmnd", config)
    failure_messages.extend(check_clusterblast_files(cluster_defs, protein_seqs, db_file, logging_only=logging_only))

    return failure_messages


def run_knownclusterblast_on_record(record: Record, options: ConfigType) -> GeneralResults:
    """ Run knownclusterblast on the given record

        Arguments:
            record: the record to analyse
            options: antismash Config

        Returns:
            an instance of GeneralResults with a result for each cluster in the record
    """
    logging.info('Running known cluster search')
    clusters, proteins = load_clusterblast_database(searchtype="knownclusterblast")
    return perform_knownclusterblast(options, record, clusters, proteins)


def perform_knownclusterblast(options: ConfigType, record: Record,
                              reference_clusters: Dict[str, ReferenceCluster],
                              proteins: Dict[str, Protein]) -> GeneralResults:
    """ Run BLAST on gene cluster proteins of each cluster, parse output and
        return result rankings for each cluster

        Only compares clusters to known clusters from the MIBiG database

        Arguments:
            options: antismash Config
            record: the Record to analyse
            clusters: a dictionary mapping reference cluster name to ReferenceCluster
            proteins: a dictionary mapping reference protein name to Protein

        Returns:
            a GeneralResults instance storing results for all clusters in the
            record
    """
    logging.debug("Running DIAMOND knowncluster searches..")
    version = find_latest_database_version(os.path.join(options.database_dir, "knownclusterblast"))
    results = GeneralResults(record.id, search_type="knownclusterblast", data_version=version)

    blastoutput = run_diamond_on_all_regions(record.get_regions(), _get_datafile_path('proteins', options))
    clusters_by_number, _ = parse_all_clusters(blastoutput, record,
                                               min_seq_coverage=40,
                                               min_perc_identity=45)

    core_gene_accessions = get_core_gene_ids(record)
    for region in record.get_regions():
        region_number = region.get_region_number()
        cluster_names_to_queries = clusters_by_number.get(region_number, {})
        ranking = score_clusterblast_output(reference_clusters,
                                            core_gene_accessions,
                                            cluster_names_to_queries)
        # store results
        region_result = RegionResult(region, ranking, proteins, "knownclusterblast")
        results.add_region_result(region_result, reference_clusters, proteins)

        write_clusterblast_output(options, record, region_result, proteins,
                                  searchtype="knownclusterblast")
    results.mibig_entries = mibig_protein_homology(blastoutput, record,
                                                   reference_clusters)
    return results


def mibig_protein_homology(blastoutput: str, record: Record, clusters: Dict[str, ReferenceCluster]
                           ) -> Dict[int, Dict[str, List[MibigEntry]]]:
    """ Constructs a mapping of clusters and genes to MiBiG hits

        Arguments:
            blastoutput: a string containing blast-formatted results
            record: the record to analyse
            clusters: a dictionary mapping reference cluster name to ReferenceCluster

        Returns:
            a dict of dicts of lists, accessed by:
                mibig_entries[cluster_number][gene_accession] = list of MibigEntry
    """
    _, queries_by_region = parse_all_clusters(blastoutput, record,
                                              min_seq_coverage=20,
                                              min_perc_identity=20)
    mibig_entries = {}

    for region in record.get_regions():
        region_number = region.get_region_number()
        queries = queries_by_region.get(region_number, {})
        cluster_entries = {}
        for cluster_protein in queries.values():
            protein_name = cluster_protein.id
            protein_entries = []
            for subject in cluster_protein.subjects.values():
                mibig_id, mibig_cluster_number = subject.genecluster.rsplit('_c', 1)
                entry = MibigEntry(subject.locus_tag, subject.annotation,
                                   mibig_id, int(mibig_cluster_number),
                                   clusters[subject.genecluster].cluster_type,
                                   subject.perc_ident, subject.blastscore,
                                   subject.perc_coverage, subject.evalue)
                protein_entries.append(entry)
            cluster_entries[protein_name] = protein_entries
        if cluster_entries:
            mibig_entries[region_number] = cluster_entries
    return mibig_entries
