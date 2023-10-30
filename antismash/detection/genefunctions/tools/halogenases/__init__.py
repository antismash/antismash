# License: GNU Affero General Public License v3 or later
# A copy of GNU AGPL v3 should have been included in this software package in LICENSE.txt.

""" Identification and classification of halogenases """

import glob
import os
from typing import Any, Optional

from antismash.config import ConfigType
from antismash.common import hmmer

from ..core import Tool
from .flavin_dependent import substrate_analysis as fdh
from .flavin_dependent.substrate_analysis import HalogenaseResults, classify

TOOL_NAME = "halogenases"


def check_prereqs(options: ConfigType) -> list[str]:
    """ Check all prerequisites are satisfied
    """
    failure_messages = []

    for binary_name in ["hmmscan", "hmmsearch", "hmmpress"]:
        if binary_name not in options.executables:
            failure_messages.append(f"Failed to locate file: {binary_name!r}")
    failure_messages.extend(prepare_data(logging_only=True))
    return failure_messages


def prepare_data(logging_only: bool = False) -> list[str]:
    """ Ensures tool data is fully prepared

        Arguments:
            logging_only: whether to return error messages instead of raising exceptions

        Returns:
            a list of error messages (only if logging_only is True)
    """
    failures = []
    aggregate_file = fdh.COMBINED_FDH_PROFILES
    individual_files = []
    for filename in glob.iglob(os.path.join(os.path.dirname(aggregate_file), "*.hmm")):
        if filename == aggregate_file:
            continue
        individual_files.append(filename)

    failures.extend(hmmer.aggregate_profiles(aggregate_file, individual_files, return_not_raise=logging_only))

    return failures


def regenerate_results(data: Optional[dict[str, Any]]) -> Optional[HalogenaseResults]:
    """ Regenerates the results for the tool from the given JSON """
    if not data:
        return None
    return HalogenaseResults.from_json(data)


TOOL = Tool(
    name=TOOL_NAME,
    check_prereqs=check_prereqs,
    classify=classify,
    prepare_data=prepare_data,
)
