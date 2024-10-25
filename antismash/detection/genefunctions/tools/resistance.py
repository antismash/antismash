# License: GNU Affero General Public License v3 or later
# A copy of GNU AGPL v3 should have been included in this software package in LICENSE.txt.

""" Finds resistance genes with the Resfam core database
    as supplied by: http://www.dantaslab.org/resfams.

    Further information about the generation of the Resfam database can be found
    in:

    Gibson MK, Forsberg KJ, Dantas G.
    Improved annotation of antibiotic resistance functions reveals microbial
        resistomes cluster by ecology.
    The ISME Journal. 2014, doi:ISMEJ.2014.106
"""

import os
from typing import Any, Iterable, Optional

from antismash.common import hmmer
from antismash.common.secmet.features import CDSFeature
from antismash.common.secmet.qualifiers.gene_functions import GeneFunction
from antismash.config import ConfigType, get_config

from .core import HMMFunctionResults as Results, Tool, scan_profiles_for_functions

TOOL_NAME = "resist"


def check_prereqs(options: ConfigType) -> list[str]:
    """ Check all prerequisites are satisfied
    """
    failure_messages = []

    for binary_name in ['hmmscan', 'hmmpress']:
        if binary_name not in options.executables:
            failure_messages.append(f"Failed to locate file: {binary_name!r}")
    failure_messages.extend(prepare_data(logging_only=True))
    return failure_messages


def classify(cds_features: Iterable[CDSFeature], options: ConfigType,
             ) -> Results:
    """ Finds any resistance genes within the given CDS features

        Arguments:
            cds_features: the CDS features to check
            options: the antiSMASH config

        Returns:
            a results instance with best hits and classification for each CDS
    """
    hmm_file = os.path.join(options.database_dir, "resfam", "Resfams.hmm")
    hits = scan_profiles_for_functions(cds_features, hmm_file, hmmscan_opts=["--cut_ga"])
    mapping = {cds_name: GeneFunction.RESISTANCE for cds_name in hits}
    return Results(tool=TOOL_NAME, best_hits=hits, function_mapping=mapping)


def prepare_data(logging_only: bool = False) -> list[str]:
    """ Ensures tool data is fully prepared

        Arguments:
            logging_only: whether to return error messages instead of raising exceptions

        Returns:
            a list of error messages (only if logging_only is True)
    """
    database = os.path.join(get_config().database_dir, 'resfam', 'Resfams.hmm')

    failures = []
    # account for database directories mounted into docker containers
    if "mounted_at_runtime" not in database:
        failures.extend(hmmer.ensure_database_pressed(database, return_not_raise=logging_only))

    return failures


def regenerate_results(data: Optional[dict[str, Any]]) -> Optional[Results]:
    """ Regenerates the tool's results from the given data """
    if not data:
        return None
    return Results.from_json(data)


TOOL = Tool(
    name=TOOL_NAME,
    check_prereqs=check_prereqs,
    classify=classify,
    prepare_data=prepare_data,
)
