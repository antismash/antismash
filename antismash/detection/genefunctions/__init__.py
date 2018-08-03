# License: GNU Affero General Public License v3 or later
# A copy of GNU AGPL v3 should have been included in this software package in LICENSE.txt.

""" A module to find and mark gene functions in records with a variety of tools
"""

import logging
import os
from typing import Any, Dict, List, Optional

from antismash.common import path, module_results
from antismash.common.secmet import Record
from antismash.config import ConfigType, get_config
from antismash.config.args import ModuleArgs

from .core import FunctionResults
from .tools import run_tools

NAME = "genefunctions"
SHORT_DESCRIPTION = "Gene function annotations"


class AllFunctionResults(module_results.DetectionResults):
    """ A collection of results for a variety of gene function detections """
    schema_version = 1

    def __init__(self, record_id: str) -> None:
        super().__init__(record_id)
        self._tools = {}  # type: Dict[str, FunctionResults]

    def add_tool_results(self, results: FunctionResults) -> None:
        """ Add results for a tool, tool name must be unique and user-friendly """
        assert isinstance(results, FunctionResults), type(results)
        if results.tool in self._tools:
            raise ValueError("Gene function results already exist for tool: %s" % results.tool)
        self._tools[results.tool] = results

    def add_to_record(self, record: Record) -> None:
        for results in self._tools.values():
            results.add_to_record(record)

    @staticmethod
    def from_json(json: Dict[str, Any], record: Record) -> "AllFunctionResults":
        if json.get("schema_version") != FunctionResults.schema_version:
            logging.debug("Schema version mismatch, discarding %s results", NAME)
            return None
        if record.id != json.get("record_id"):
            logging.debug("Record ID mismatch, discarding %s results", NAME)
            return None

        all_results = AllFunctionResults(record.id)
        for json_result in json["tools"]:
            all_results.add_tool_results(FunctionResults.from_json(json_result, record))

        return all_results

    def to_json(self) -> Dict[str, Any]:
        return {"schema_version": self.schema_version,
                "record_id": self.record_id,
                "tools": [res.to_json() for res in self._tools.values()]
                }


def get_arguments() -> ModuleArgs:
    """ Build and return arguments. No extra options beyond a switch to enable """
    args = ModuleArgs('Additional analysis', 'genefunctions', enabled_by_default=True)
    return args


def check_options(_options: ConfigType) -> List[str]:
    """ Checks options for conflicts.
        No extra options, so they can't have conflicts.
    """
    return []


def check_prereqs() -> List[str]:
    """Check for prerequisites
    """
    failure_messages = []
    database_dir = os.path.join(get_config().database_dir, 'resfam', 'Resfams.hmm')
    if path.locate_file(database_dir) is None:
        failure_messages.append('Failed to locate Resfam database in %s' % database_dir)
    return failure_messages


def is_enabled(options: ConfigType) -> bool:
    """ Should the module be run with these options """
    return not options.minimal or options.genefunctions_enabled


def regenerate_previous_results(previous: Dict[str, Any], record: Record,
                                _options: ConfigType) -> Optional[AllFunctionResults]:
    """ Regenerate the previous results from JSON format. """
    if not previous:
        return None
    return AllFunctionResults.from_json(previous, record)


def run_on_record(record: Record, results: AllFunctionResults, options: ConfigType) -> AllFunctionResults:
    """ Run the detection, unless the previous results apply to the given record """
    if isinstance(results, AllFunctionResults) and results.record_id == record.id:
        return results
    results = AllFunctionResults(record.id)
    for result in run_tools(record, options):
        results.add_tool_results(result)
    return results
