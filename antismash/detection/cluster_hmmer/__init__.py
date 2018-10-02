# License: GNU Affero General Public License v3 or later
# A copy of GNU AGPL v3 should have been included in this software package in LICENSE.txt.

""" PFAM anotation for only clusters """

import logging
import os
from typing import Any, Dict, List, Optional

from antismash.common import hmmer, path, pfamdb
from antismash.common.secmet import Record
from antismash.config import get_config, ConfigType
from antismash.config.args import ModuleArgs

NAME = "cluster_hmmer"
SHORT_DESCRIPTION = "Cluster PFAM anotation"

MIN_SCORE = 0.
MAX_EVALUE = 0.01


def get_arguments() -> ModuleArgs:
    """ Builds the module args """
    args = ModuleArgs('Cluster HMMer options', 'clusterhmmer')
    args.add_analysis_toggle('clusterhmmer',
                             dest='clusterhmmer',
                             action='store_true',
                             default=False,
                             help="Run a cluster-limited HMMer analysis.")
    args.add_option('pfamdb-version',
                    dest='pfamdb_version',
                    type=str,
                    default='latest',
                    help="PFAM database version number (e.g. 27.0) (default: %(default)s).")
    return args


def is_enabled(options: ConfigType) -> bool:
    """  Uses the supplied options to determine if the module should be run """
    return options.clusterhmmer


def check_prereqs() -> List[str]:
    """ Ensure at least one database exists and is valid """
    failure_messages = []
    for binary_name in ['hmmscan']:
        if not path.locate_executable(binary_name):
            failure_messages.append("Failed to locate executable: %r" % binary_name)

    data_dir = get_config().database_dir
    try:
        version = pfamdb.find_latest_database_version(data_dir)
    except ValueError as err:
        failure_messages.append(str(err))
        return failure_messages

    data_path = os.path.join(data_dir, "pfam", version)
    failure_messages.extend(pfamdb.check_db(data_path))
    return failure_messages


def check_options(options: ConfigType) -> List[str]:
    """ Check the requested PFAM database exists """
    database_version = options.clusterhmmer_pfamdb_version
    pfam_dir = os.path.join(options.database_dir, "pfam")
    if database_version == "latest":
        database_version = pfamdb.find_latest_database_version(options.database_dir)
    return pfamdb.check_db(os.path.join(pfam_dir, database_version))


def regenerate_previous_results(previous: Dict[str, Any], record: Record,
                                _options: ConfigType) -> Optional[hmmer.HmmerResults]:
    """ Rebuild previous results """
    if not previous:
        return None
    return hmmer.HmmerResults.from_json(previous, record, MAX_EVALUE, MIN_SCORE)


def run_on_record(record: Record, results: Optional[hmmer.HmmerResults],
                  options: ConfigType) -> hmmer.HmmerResults:
    """ Run hmmsearch against PFAM for all CDS features within the record """

    if options.clusterhmmer_pfamdb_version == "latest":
        database_version = pfamdb.find_latest_database_version(options.database_dir)
    else:
        database_version = options.clusterhmmer_pfamdb_version

    if results:
        previous_db = pfamdb.get_db_version_from_path(results.database)
        # same version requested, so reuse the results
        if database_version == previous_db:
            return results
        else:
            logging.debug("Replacing clusterhmmer results from %s with %s",
                          previous_db, database_version)

    logging.info('Running cluster PFAM search')

    features = []
    for region in record.get_regions():
        features.extend(list(region.cds_children))
    database = os.path.join(options.database_dir, 'pfam', database_version, 'Pfam-A.hmm')
    return hmmer.run_hmmer(record, features, MAX_EVALUE, MIN_SCORE, database, "clusterhmmer")
