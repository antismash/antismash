# License: GNU Affero General Public License v3 or later
# A copy of GNU AGPL v3 should have been included in this software package in LICENSE.txt.

""" Terpene analysis module
"""

import logging
from typing import Any, Optional

from antismash.common import hmmer, path
from antismash.common.secmet import Record
from antismash.config import ConfigType
from antismash.config.args import ModuleArgs

from .data_loader import load_hmm_properties
from .results import TerpeneResults
from .terpene_analysis import analyse_cluster
from .html_output import generate_html, will_handle

NAME = "terpene"
SHORT_DESCRIPTION = "Terpene analysis"


def get_arguments() -> ModuleArgs:
    """ Return the args for the terpene module """
    args = ModuleArgs('Advanced options', 'terpene', enabled_by_default=True)
    return args


def is_enabled(options: ConfigType) -> bool:
    """ Whether the module is enabled """
    return not options.minimal or options.terpene_enabled


def check_options(_options: ConfigType) -> list[str]:
    """ No options to check at the moment """
    return []


def check_prereqs(options: ConfigType) -> list[str]:
    """ Check the prerequisites.
            hmmscan: domain detection
            hmmpress: compressing the hmm files

        Returns:
            a list of strings describing any errors, if they occurred
    """
    failure_messages = []
    for binary_name in ["hmmsearch", "hmmpress"]:
        if binary_name not in options.executables:
            failure_messages.append(f"Failed to locate executable for {binary_name!r}")

    # no point checking the data if we can't use it
    if failure_messages:
        return failure_messages

    failure_messages.extend(prepare_data(logging_only=True))

    return failure_messages


def prepare_data(logging_only: bool = False) -> list[str]:
    """ Ensures packaged data is fully prepared

        Arguments:
            logging_only: whether to return error messages instead of raising exceptions

        Returns:
            a list of error messages (only if logging_only is True)
    """
    failure_messages = []

    # Check if hmm_properties.json and compound_keys.json are well-formatted
    try:
        load_hmm_properties()
    except ValueError as err:
        if not logging_only:
            raise
        failure_messages.append(str(err))

    # the hmm files that need to be present in data
    hmm_names = ["main", "PT_FPPS_like", "PT_phytoene_like", "T1TS", "T2TS"]
    hmm_files = [path.get_full_path(__file__, "data", name + "_profiles.hmm") for name in hmm_names]

    # the path to the combined data file of all hmms
    all_hmms = path.get_full_path(__file__, 'data', 'all_profiles.hmm')
    failure_messages.extend(hmmer.aggregate_profiles(all_hmms, hmm_files, return_not_raise=logging_only))

    return failure_messages


def regenerate_previous_results(previous: dict[str, Any], record: Record,
                                _options: ConfigType) -> Optional[TerpeneResults]:
    """ Regenerate the previous results from JSON format.

        Arguments:
            previous: the previous results as from JSON
            record: the Record these previous results were originally created from
            options: the current antismash config object

        Returns:
            an instance of the module's ModuleResults implementation,
            or None if the current options require the analysis to be rerun or cannot be regenerated
    """
    # if there isn't anything to work with, just return None
    if not previous:
        return None
    return TerpeneResults.from_json(previous, record)


def run_on_record(record: Record, results: TerpeneResults, _options: ConfigType) -> TerpeneResults:
    """ Run the analysis, unless the previous results apply to the given record

        Arguments:
            record: the Record being analysed
            results: an existing instance of the module's ModuleResults implementation (or None)
            options: the current antismash config object

        Returns:
            a TerpeneResults instance
    """
    if isinstance(results, TerpeneResults) and results.record_id == record.id:
        return results

    results = TerpeneResults(record.id)

    terpene_clusters = [cluster for cluster in record.get_protoclusters()
                        if cluster.product_category == 'terpene']
    if not terpene_clusters:
        logging.debug("No terpene clusters to analyse")
        return results

    logging.info("Analysing terpene clusters")
    for cluster in terpene_clusters:
        results.cluster_predictions[cluster.get_protocluster_number()] = analyse_cluster(cluster)

    return results
