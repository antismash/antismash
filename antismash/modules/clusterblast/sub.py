# License: GNU Affero General Public License v3 or later
# A copy of GNU AGPL v3 should have been included in this software package in LICENSE.txt.

""" Compares subsets of clusters to a reference data set of other subclusters
"""

import functools
import logging
from typing import Dict, List

from helperlibs.wrappers.io import TemporaryDirectory

from antismash.common import path, subprocessing
from antismash.common.secmet import Record
from antismash.config import ConfigType

from .core import write_raw_clusterblastoutput, write_fastas_with_all_genes, \
                  blastparse, score_clusterblast_output, \
                  run_blast, load_clusterblast_database, get_core_gene_ids
from .results import RegionResult, GeneralResults
from .data_structures import ReferenceCluster, Protein


def _get_datafile_path(filename: str) -> str:
    """ A simple helper to get the full path of subclusterblast datafile """
    return path.get_full_path(__file__, 'data', 'sub', filename)


def check_sub_prereqs(_options: ConfigType) -> List[str]:
    """ Check if all required applications and datafiles are present.
        options is irrelevant here
    """
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


def run_subclusterblast_on_record(record: Record, options: ConfigType) -> GeneralResults:
    """ Loads reference databases and performs subclusterblast analysis

        Arguments:
            record: the Record to analyse
            options: antismash Config

        Returns:
            a GeneralResults with results for each cluster in the record
    """
    logging.info('Running subcluster search')
    clusters, proteins = load_clusterblast_database(record, "subclusterblast")
    return perform_subclusterblast(options, record, clusters, proteins)


def _run_blast_helper(database: str, index: int) -> str:
    """ A simple wrapper of run_blast() to reverse arg order to allow for
        functools.partial to simplify run_clusterblast_processes()

        Cannot be a locally defined function because Pool requires it to be
        pickable.

        Arguments:
            database: the database to search in
            index: the index of the current process, used to generate the
                    filename to use

        Returns:
            the name of the output file created by run_blast()
    """
    return run_blast("input%d.fasta" % index, database)


def run_clusterblast_processes(options: ConfigType) -> None:
    """ Run blast in parallel, creates `options.cpu` number of files with
        the name format "inputN.fasta"

        Arguments:
            options: antismash Config

        Returns:
            None
    """
    database = _get_datafile_path('subclusterprots.fasta')
    # set the first arg to always be database
    partial = functools.partial(_run_blast_helper, database)
    # run in parallel
    subprocessing.parallel_function(partial, [[i] for i in range(options.cpus)],
                                    cpus=options.cpus)


def read_clusterblast_output(options: ConfigType) -> str:
    """ Builds a single output from the results from the distributed blast run

        Arguments:
            options: antismash Config

        Returns:
            a string containing all blast run output
    """
    blastoutput = []
    for i in range(options.cpus):
        with open("input%d.out" % i, "r") as handle:
            output = handle.read()
        blastoutput.append(output)
    return "".join(blastoutput)


def perform_subclusterblast(options: ConfigType, record: Record, clusters: Dict[str, ReferenceCluster],
                            proteins: Dict[str, Protein]) -> GeneralResults:
    """ Run BLAST on gene cluster proteins of each cluster, parse output and
        return result rankings for each cluster

        Arguments:
            options: antismash Config
            record: the Record to analyse
            clusters: a dictionary mapping reference cluster name to ReferenceCluster
            proteins: a dictionary mapping reference protein name to Protein

        Returns:
            a GeneralResults instance storing results for all clusters in the
            record
    """
    results = GeneralResults(record.id, search_type="subclusterblast")
    with TemporaryDirectory(change=True):
        allcoregenes = get_core_gene_ids(record)
        for region in record.get_regions():
            # prepare and run diamond
            write_fastas_with_all_genes([region], "input.fasta",
                                        partitions=options.cpus)
            run_clusterblast_processes(options)
            blastoutput = read_clusterblast_output(options)
            write_raw_clusterblastoutput(options.output_dir, blastoutput, prefix="subclusterblast")
            # parse and score diamond results
            _, cluster_names_to_queries = blastparse(blastoutput, record,
                                                     min_seq_coverage=40,
                                                     min_perc_identity=45)
            ranking = score_clusterblast_output(clusters, allcoregenes, cluster_names_to_queries)
            logging.debug("Cluster at %s has %d subclusterblast results", region.location, len(ranking))
            # store results
            region_result = RegionResult(region, ranking, proteins, "subclusterblast")
            results.add_region_result(region_result, clusters, proteins)
    return results
