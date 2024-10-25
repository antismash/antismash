# License: GNU Affero General Public License v3 or later
# A copy of GNU AGPL v3 should have been included in this software package in LICENSE.txt.

""" The classification section of the smCOG module. Categorises gene function
    according to a curated set of HMM profiles.
"""
from typing import Any, Iterable, Optional

from antismash.common import hmmer, path
from antismash.common.secmet import CDSFeature, GeneFunction
from antismash.config import ConfigType

from .core import HMMFunctionResults as Results, HMMHit, Tool, scan_profiles_for_functions

DATABASE = path.get_full_path(__file__, "data", "smcogs.hmm")
TOOL_NAME = "smcogs"


def build_function_mapping() -> dict[str, GeneFunction]:
    """ Load the smCOG gene function mapping from a file.

        Returns:
             a dictionary mapping smCOG id to its matching gene function
    """
    mapping = {
        'B': GeneFunction.ADDITIONAL,  # 'biosynthetic-additional',
        'T': GeneFunction.TRANSPORT,  # 'transport',
        'R': GeneFunction.REGULATORY,  # 'regulatory',
    }
    annotations = {}
    with open(path.get_full_path(__file__, "data", "cog_annotations.txt"), "r", encoding="utf-8") as handle:
        lines = handle.readlines()
    for line in lines:
        cog, _desc, key = line.strip().split('\t', 3)
        annotations[cog] = mapping.get(key, GeneFunction.OTHER)

    return annotations


def check_prereqs(options: ConfigType) -> list[str]:
    """ Check all prerequisites are satisfied
    """
    failure_messages = []

    for binary_name in ['hmmscan', 'hmmpress']:
        if binary_name not in options.executables:
            failure_messages.append(f"Failed to locate file: {binary_name!r}")
    failure_messages.extend(prepare_data(logging_only=True))
    return failure_messages


def classify(cds_features: Iterable[CDSFeature],
             options: ConfigType) -> Results:  # pylint: disable=unused-argument  # classify() is an API
    """ Finds possible classifications for the provided CDS features.

        Arguments:
            cds_features: a list of CDSFeatures to classify
            options: the current run options

        Returns:
            a results instance with best hits and classification for each CDS
    """

    hits: dict[str, HMMHit] = scan_profiles_for_functions(cds_features, DATABASE, hmmscan_opts=["-E", "1E-16"])
    ids_to_function = build_function_mapping()
    cds_name_to_function = {}
    for cds_name, hit in hits.items():
        # pull out the identifier by itself
        smcog_id, description = hit.reference_id.split(":", 1)
        # remove extraneous info
        hit.reference_id = smcog_id
        hit.description = description.replace("_", " ")
        cds_name_to_function[cds_name] = ids_to_function[smcog_id]

    return Results(tool=TOOL_NAME, best_hits=hits, function_mapping=cds_name_to_function)


def prepare_data(logging_only: bool = False) -> list[str]:
    """ Ensures tool data is fully prepared

        Arguments:
            logging_only: whether to return error messages instead of raising exceptions

        Returns:
            a list of error messages (only if logging_only is True)
    """
    failures = []
    failures.extend(hmmer.ensure_database_pressed(DATABASE, return_not_raise=logging_only))
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
