# License: GNU Affero General Public License v3 or later
# A copy of GNU AGPL v3 should have been included in this software package in LICENSE.txt.

""" The classification section of the smCOG module. Categorises gene function
    according to a curated set of HMM profiles.
"""

from dataclasses import dataclass
from typing import Any, Iterable, Optional, Self

from antismash.common import hmmer, json, path
from antismash.common.secmet import CDSFeature, ECGroup, GeneFunction
from antismash.config import ConfigType

from .core import HMMFunctionResults as Results, HMMHit, Tool, scan_profiles_for_functions

DATABASE = path.get_full_path(__file__, "data", "smcogs.hmm")
METADATA = path.get_full_path(__file__, "data", "cog_annotations.json")
TOOL_NAME = "smcogs"


@dataclass
class SmcogProfile:
    """ A container for all information regarding a profile """
    id: str
    description: str
    function: GeneFunction
    groups: list[ECGroup]
    subfunctions: list[str]

    @classmethod
    def from_json(cls, data: dict[str, Any]) -> Self:
        """ Reconstructs an instance from JSON """
        return cls(
            id=data["id"],
            description=data["description"],
            function=GeneFunction(data["function"]),
            groups=[ECGroup(g) for g in data["groups"]],
            subfunctions=data["subfunctions"],
        )


@dataclass
class SmcogMetadata:
    """ A container for SMCoG data """
    version: str
    profiles: list[SmcogProfile]

    @classmethod
    def from_json(cls, data: dict[str, Any]) -> Self:
        """ Reconstructs an instance from JSON """
        return cls(
            version=data["version"],
            profiles=[SmcogProfile.from_json(p) for p in data["profiles"]],
        )


def _load_profiles() -> dict[str, SmcogProfile]:
    """ Load the smCOG metadata from a file.

        Returns:
             a dictionary mapping smCOG id to its matching SmcogProfile
    """
    mapping = {}
    with open(METADATA, "rb") as handle:
        data = json.loads(handle.read())
    metadata = SmcogMetadata.from_json(data)

    for profile in metadata.profiles:
        mapping[profile.id] = profile

    return mapping


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
    profiles = _load_profiles()
    function_mapping = {}
    group_mapping = {}
    subfunction_mapping = {}
    for cds_name, hit in hits.items():
        # pull out the identifier by itself
        smcog_id = hit.reference_id.split(":", 1)[0]
        hit.reference_id = smcog_id
        profile = profiles[smcog_id]
        hit.description = profile.description
        function_mapping[cds_name] = profile.function
        group_mapping[cds_name] = profile.groups
        subfunction_mapping[cds_name] = profile.subfunctions

    return Results(tool=TOOL_NAME, best_hits=hits, function_mapping=function_mapping,
                   group_mapping=group_mapping, subfunction_mapping=subfunction_mapping)


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
