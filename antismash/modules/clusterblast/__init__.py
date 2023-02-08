# License: GNU Affero General Public License v3 or later
# A copy of GNU AGPL v3 should have been included in this software package in LICENSE.txt.

"""ClusterBlast comparative gene cluster analysis"""

import logging
import os
from typing import Any, Dict, List, Optional

from antismash.common.secmet import Record
from antismash.config import get_config, ConfigType
from antismash.config.args import ModuleArgs

from .core import (
    check_clusterblast_files,
    internal_homology_blast,
    load_clusterblast_database,
)
from .clusterblast import perform_clusterblast
from .html_output import generate_html, will_handle
from .known import run_knownclusterblast_on_record, check_known_prereqs, prepare_known_data
from .results import ClusterBlastResults, get_result_limit
from .sub import run_subclusterblast_on_record, check_sub_prereqs

NAME = "clusterblast"
SHORT_DESCRIPTION = "comparative gene cluster analysis"


def get_arguments() -> ModuleArgs:
    """ Builds the args for the clusterblast module """
    args = ModuleArgs('ClusterBlast options', 'cb')
    args.add_analysis_toggle('general',
                             dest='general',
                             action='store_true',
                             default=False,
                             help="Compare identified clusters against a "
                                  "database of antiSMASH-predicted clusters.")
    args.add_analysis_toggle('subclusters',
                             dest='subclusters',
                             action='store_true',
                             default=False,
                             help="Compare identified clusters against known "
                                  "subclusters responsible for synthesising "
                                  "precursors.")
    args.add_analysis_toggle('knownclusters',
                             dest='knownclusters',
                             action='store_true',
                             default=False,
                             help="Compare identified clusters against known "
                                  "gene clusters from the MIBiG database.")
    args.add_option('nclusters',
                    dest='nclusters',
                    metavar="count",
                    type=int,
                    default=10,
                    help="Number of clusters from ClusterBlast to display,"
                         f" cannot be greater than {get_result_limit()}. (default: %(default)s)")
    args.add_option('min-homology-scale',
                    dest='min_homology_scale',
                    metavar="LIMIT",
                    type=float,
                    default=0.0,
                    help="A minimum scaling factor for the query BGC in ClusterBlast results."
                         " Valid range: 0.0 - 1.0.   "
                         " Warning: some homologous genes may no longer be visible!"
                         " (default: %(default)s)")
    return args


def is_enabled(options: ConfigType) -> bool:
    """  Uses the supplied options to determine if the module should be run
    """
    return options.cb_general or options.cb_knownclusters or options.cb_subclusters


def check_options(options: ConfigType) -> List[str]:
    """ Checks that extra options are valid """
    if options.cb_nclusters > get_result_limit():
        return [f"nclusters of {options.cb_nclusters} is over limit of {get_result_limit()}"]
    return []


def regenerate_previous_results(previous: Dict[str, Any], record: Record,
                                _options: ConfigType) -> Optional[ClusterBlastResults]:
    """ Regenerates previous results """
    if not previous:
        logging.debug("No previous clusterblast results to reuse")
        return None
    return ClusterBlastResults.from_json(previous, record)


def check_prereqs(options: ConfigType) -> List[str]:
    "Check if all required applications are around"
    _required_binaries = [
        'blastp',
        'makeblastdb',
        'diamond'
    ]

    failure_messages = []
    for binary_name in _required_binaries:
        if binary_name not in options.executables:
            failure_messages.append(f"Failed to locate file: {binary_name!r}")

    if "diamond" not in get_config().executables:
        failure_messages.append("cannot check clusterblast databases, no diamond executable present")
        return failure_messages

    failure_messages.extend(prepare_data(logging_only=True))

    failure_messages.extend(check_known_prereqs(options))
    failure_messages.extend(check_sub_prereqs(options))

    return failure_messages


def prepare_data(logging_only: bool = False) -> List[str]:
    """ Prepare the databases. """
    failure_messages = []
    # known
    failure_messages.extend(prepare_known_data(logging_only))

    # general
    clusterblastdir = os.path.join(get_config().database_dir, "clusterblast")
    if "mounted_at_runtime" in clusterblastdir:  # can't prepare these
        return failure_messages
    cluster_defs = os.path.join(clusterblastdir, 'clusters.txt')
    protein_seqs = os.path.join(clusterblastdir, "proteins.fasta")
    db_file = os.path.join(clusterblastdir, "proteins.dmnd")

    # check the DBv3 region info exists instead of single cluster numbers
    with open(protein_seqs, encoding="utf-8") as handle:
        sample = handle.readline()
    if "-" not in sample.split("|", 3)[1]:
        failure_messages.append("clusterblast database out of date, update with download-databases")
        # and don't bother pressing them
        return failure_messages

    failure_messages.extend(check_clusterblast_files(cluster_defs, protein_seqs, db_file, logging_only=logging_only))

    return failure_messages


def run_on_record(record: Record, results: Optional[ClusterBlastResults],
                  options: ConfigType) -> ClusterBlastResults:
    """ Runs the specified clusterblast variants over the record """
    if not results:
        results = ClusterBlastResults(record.id)
        results.internal_homology_groups = internal_homology_blast(record)
    if options.cb_general and not results.general:
        logging.info('Running ClusterBlast')
        clusters, proteins = load_clusterblast_database()
        results.general = perform_clusterblast(options, record, clusters, proteins)
    if options.cb_subclusters and not results.subcluster:
        results.subcluster = run_subclusterblast_on_record(record, options)
    if options.cb_knownclusters and not results.knowncluster:
        results.knowncluster = run_knownclusterblast_on_record(record, options)
    return results
