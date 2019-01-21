# License: GNU Affero General Public License v3 or later
# A copy of GNU AGPL v3 should have been included in this software package in LICENSE.txt.

"""ClusterBlast comparative gene cluster analysis"""

import logging
import os
from typing import Any, Dict, List, Optional

from antismash.common import path
from antismash.common.secmet import Record
from antismash.config import get_config, ConfigType
from antismash.config.args import ModuleArgs

from .core import load_clusterblast_database, internal_homology_blast
from .clusterblast import perform_clusterblast
from .html_output import generate_html, will_handle
from .known import run_knownclusterblast_on_record, check_known_prereqs
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
                         " cannot be greater than %d. (default: %%(default)s)" % get_result_limit())
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
        return ["nclusters of %d is over limit of %d" % (
                    options.cb_nclusters, get_result_limit())]
    return []


def regenerate_previous_results(previous: Dict[str, Any], record: Record,
                                _options: ConfigType) -> Optional[ClusterBlastResults]:
    """ Regenerates previous results """
    if not previous:
        logging.debug("No previous clusterblast results to reuse")
        return None
    return ClusterBlastResults.from_json(previous, record)


def check_prereqs() -> List[str]:
    "Check if all required applications are around"
    options = get_config()
    # Tuple is ( binary_name, optional)
    _required_binaries = [
        ('blastp', False),
        ('makeblastdb', False),
        ('diamond', False),
    ]

    _required_files = [
        ('geneclusterprots.dmnd', False),
        ('geneclusterprots.fasta', False),
        ('geneclusters.txt', False),
    ]

    clusterblastdir = os.path.join(options.database_dir, "clusterblast")

    failure_messages = []
    for binary_name, optional in _required_binaries:
        if path.locate_executable(binary_name) is None and not optional:
            failure_messages.append("Failed to locate file: %r" % binary_name)

    for file_name, optional in _required_files:
        if path.locate_file(os.path.join(clusterblastdir, file_name)) is None and not optional:
            failure_messages.append("Failed to locate file: %r" % file_name)

    failure_messages.extend(check_known_prereqs(options))
    failure_messages.extend(check_sub_prereqs(options))
    return failure_messages


def run_on_record(record: Record, results: Optional[ClusterBlastResults],
                  options: ConfigType) -> ClusterBlastResults:
    """ Runs the specified clusterblast variants over the record """
    if not results:
        results = ClusterBlastResults(record.id)
        results.internal_homology_groups = internal_homology_blast(record)
    if options.cb_general and not results.general:
        logging.info('Running ClusterBlast')
        clusters, proteins = load_clusterblast_database(record)
        results.general = perform_clusterblast(options, record, clusters, proteins)
    if options.cb_subclusters and not results.subcluster:
        results.subcluster = run_subclusterblast_on_record(record, options)
    if options.cb_knownclusters and not results.knowncluster:
        results.knowncluster = run_knownclusterblast_on_record(record, options)
    return results
