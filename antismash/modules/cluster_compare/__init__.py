# License: GNU Affero General Public License v3 or later
# A copy of GNU AGPL v3 should have been included in this software package in LICENSE.txt.

""" A module for comparing clusters """

import json
import os
from typing import Any, Dict, List, Optional

import jsonschema

from antismash.common import path
from antismash.config import ConfigType, get_config
from antismash.config.args import ModuleArgs, MultipleFullPathAction
from antismash.common.secmet import Record
from antismash.common.subprocessing.diamond import check_diamond_files

from .analysis import run
from .data_structures import DBConfig
from .html_output import generate_html, will_handle, generate_javascript_data
from .results import ClusterCompareResults

NAME = "cluster-compare"
SHORT_DESCRIPTION = "cluster comparison"

_MIBIG_CONFIG = path.get_full_path(__file__, "data", "mibig.json")


def get_arguments() -> ModuleArgs:
    """ Constucts module arguments
    """
    args = ModuleArgs("ClusterCompare options", "cc")
    args.add_analysis_toggle("mibig",
                             dest='mibig',
                             action='store_true',
                             default=False,
                             help="Run a comparison against the MIBiG dataset")
    args.add_option("custom-dbs",
                    dest="custom_dbs",
                    metavar="FILE1,FILE2,...",
                    action=MultipleFullPathAction,
                    default=[],
                    help="A comma separated list of database config files to run with")
    return args


def _get_all_databases(options: ConfigType, defaults: bool = False) -> List[DBConfig]:
    """ Returns a list of DBConfig instances for each database specified in options

        Does no validation, so should only be used after validation steps.

        Arguments:
            options: the current config
            defaults: whether to ignore the state of built-in database options
                      and prepare them anyway
        Returns:
            a list of DBConfig instances
    """
    dbs = []
    if options.cc_mibig or defaults:
        dbs.append(DBConfig.from_file(_MIBIG_CONFIG, options.database_dir))
    for custom_db in options.cc_custom_dbs:
        dbs.append(DBConfig.from_file(custom_db, options.database_dir))
    return dbs


def check_options(options: ConfigType) -> List[str]:
    """ Check the options of this module for any conflicting or invalid values.

        Arguments:
            options: the options parsed by the main entry point as an
                     antismash.Config object

        Returns:
            a list of strings describing any errors, if they exist
    """
    errors = []
    with open(path.get_full_path(__file__, "data", "schema.json"), encoding="utf-8") as handle:
        schema = json.load(handle)
    for db in options.cc_custom_dbs:
        try:
            with open(db, encoding="utf-8") as handle:
                setup = json.load(handle)
                try:
                    jsonschema.validate(instance=setup, schema=schema)
                except jsonschema.ValidationError as err:
                    errors.append(f"custom clustercompare database {db} is invalid: {err.message}")
                    continue
                config = DBConfig.from_json(setup, options.database_dir)
                valid = True
                for filename in ["data.json", "proteins.fasta"]:
                    full = os.path.join(config.path, filename)
                    if not path.locate_file(full):
                        errors.append(f"custom clustercompare database '{setup['name']}' missing data file {full}")
                        valid = False
                if not valid:
                    continue
        except FileNotFoundError:
            errors.append(f"could not open custom clustercompare database file {db}")
    return errors


def prepare_data(logging_only: bool = False) -> List[str]:
    """ Ensures packaged data is fully prepared

        Arguments:
            logging_only: whether to return error messages instead of raising exceptions

        Returns:
            a list of error messages (only if logging_only is True)
    """
    failure_messages = []
    options = get_config()

    for db in _get_all_databases(options, defaults=True):
        # account for database directories mounted into docker containers
        if "mounted_at_runtime" in db.path:
            continue
        cluster_defs = os.path.join(db.path, 'data.json')
        protein_seqs = os.path.join(db.path, "proteins.fasta")
        db_file = os.path.join(db.path, "proteins.dmnd")
        failure_messages.extend(check_diamond_files(cluster_defs, protein_seqs,
                                                    db_file, logging_only=logging_only))
    return failure_messages


def check_prereqs(options: ConfigType) -> List[str]:
    "Check if all required applications are around"

    failure_messages = []
    for binary_name in ["diamond"]:
        if binary_name not in options.executables:
            failure_messages.append(f"Failed to locate file: {binary_name!r}")

    if "diamond" not in get_config().executables:
        failure_messages.append("cannot check cluster_compare databases, no diamond executable present")
        return failure_messages

    failure_messages.extend(prepare_data(logging_only=True))
    return failure_messages


def is_enabled(options: ConfigType) -> bool:
    """ Returns True if the module is enabled by the given options """
    return options.cc_mibig or options.cc_custom_dbs


def regenerate_previous_results(previous: Dict[str, Any], record: Record,
                                _options: ConfigType) -> Optional[ClusterCompareResults]:
    """ Rebuild the previous run results from a JSON object into this module's
        python results class. If the current options are incompatible with the
        previous results, None should be returned.

        Arguments:
            previous: the previous results as a dictionary
            record: the Record that was used to generate the previous results
            options: an antismash.Config object
    """
    return ClusterCompareResults.from_json(previous, record)


def run_on_record(record: Record, results: Optional[ClusterCompareResults],
                  options: ConfigType) -> ClusterCompareResults:
    """ Run this module's analysis on the requested databases. If previous results
        exist for a particular database, they will be used instead of rerunning.
        If no previous results exist for a database, they will be added to any
        existing results.

        Arguments:
            record: the Record instance to analyse
            results: the previous results as generated by regenerate_previous_results()
            options: an antismash.Config object

        Returns:
            a ClusterCompareResults instance
    """
    if results is None:
        results = ClusterCompareResults(record.id, {})

    for db in _get_all_databases(options):
        # if reusing results, skip results for databases that have already been run
        if db.name in results.by_database:
            continue
        results.by_database[db.name] = run(record, db)

    return results
