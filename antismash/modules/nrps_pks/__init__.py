# License: GNU Affero General Public License v3 or later
# A copy of GNU AGPL v3 should have been included in this software package in LICENSE.txt.

""" NRPS/PKS analysis module
"""

from typing import Any, Dict, List, Optional

from antismash.common import path
from antismash.common.brawn import CacheError, ensure_alignment_cached
from antismash.common.secmet import Record
from antismash.config import ConfigType
from antismash.config.args import ModuleArgs

from .html_output import generate_html, will_handle
from .results import NRPS_PKS_Results
from .specific_analysis import specific_analysis
from .at_analysis import prepare_data as at_prepare_data
from .kr_analysis import prepare_data as kr_prepare_data
from .minowa import prepare_data as minowa_prepare_data
from .nrpys import check_prereqs as nrpys_check_prereqs
from .orderfinder import C_TERMINAL_PATH, N_TERMINAL_PATH

NAME = "nrps_pks"

SHORT_DESCRIPTION = "NRPS/PKS analysis"


def prepare_data(logging_only: bool = False) -> List[str]:
    """ Ensures packaged data is fully prepared

        Arguments:
            logging_only: whether to return error messages instead of raising exceptions

        Returns:
            a list of error messages (only if logging_only is True)
    """
    failures: List[str] = []
    # pre-cache alignments need to be cached and up to date
    expected_alignments = [
        C_TERMINAL_PATH,
        N_TERMINAL_PATH,
    ]
    cache_dir = path.get_full_path(__file__, "data")
    for fasta in expected_alignments:
        try:
            ensure_alignment_cached(fasta, cache_dir)
        except CacheError as err:
            if not logging_only:
                raise
            failures.append(str(err))
    for func in [at_prepare_data, kr_prepare_data, minowa_prepare_data]:
        failures.extend(func(logging_only=logging_only))
    return failures


def check_prereqs(options: ConfigType) -> List[str]:
    """ Check the prerequisites.
            hmmsearch: minowa
    """
    failure_messages = []
    for binary_name in ["hmmsearch"]:
        if binary_name not in options.executables:
            failure_messages.append(f"Failed to locate executable for {binary_name!r}")
    failure_messages.extend(prepare_data(logging_only=True))
    failure_messages.extend(nrpys_check_prereqs(options))
    return failure_messages


def get_arguments() -> ModuleArgs:
    """ Construct module arguments and for sub-modules """
    args = ModuleArgs('NRPS/PKS options', 'np', enabled_by_default=True)
    return args


def check_options(_options: ConfigType) -> List[str]:
    """ No options to check at the moment """
    return []


def regenerate_previous_results(json: Dict[str, Any], record: Record,
                                _options: ConfigType) -> Optional[NRPS_PKS_Results]:
    """ Regenerate results from previous run """
    return NRPS_PKS_Results.from_json(json, record)


def is_enabled(options: ConfigType) -> bool:
    """ Whether the module is enabled """
    return not options.minimal or options.nrps_pks_enabled


def run_on_record(record: Record, results: Optional[NRPS_PKS_Results], options: ConfigType) -> NRPS_PKS_Results:
    """ Analyse a record """
    if not results:
        results = NRPS_PKS_Results(record.id)
        specific_analysis(record, results, options)
    return results
