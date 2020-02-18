# License: GNU Affero General Public License v3 or later
# A copy of GNU AGPL v3 should have been included in this software package in LICENSE.txt.

"""Map Pfam domains to Gene Ontology terms"""

from typing import Any, Dict, List, Optional

from antismash.common import path
from antismash.common.secmet import Record
from antismash.config import ConfigType
from antismash.config.args import ModuleArgs
from .pfam2go import DATA_FILE, get_gos_for_pfams, Pfam2GoResults

NAME = "pfam2go"
SHORT_DESCRIPTION = "Pfam domain to Gene Ontology mapping"


def get_arguments() -> ModuleArgs:
    """ Build and return arguments. No extra options beyond a switch to enable """
    args = ModuleArgs('Additional analysis', 'pfam2go')
    args.add_analysis_toggle('--pfam2go',
                             dest='pfam2go',
                             action='store_true',
                             default=False,
                             help="Run Pfam to Gene Ontology mapping module.")
    return args


def check_options(_options: ConfigType) -> List[str]:
    """ Checks options for conflicts.
        No extra options, so they can't have conflicts.
    """
    return []


def check_prereqs(_options: ConfigType) -> List[str]:
    """Check for prerequisites
        data file: mapping file for Pfam to Gene Ontology mapping
    """
    failure_messages = []
    if path.locate_file(DATA_FILE) is None:
        failure_messages.append('Failed to locate Pfam to Gene Ontology mapping file')
    return failure_messages


def is_enabled(options: ConfigType) -> bool:
    """ Should the module be run with these options """
    return options.pfam2go


def regenerate_previous_results(previous: Dict[str, Any], record: Record,
                                _options: ConfigType) -> Optional[Pfam2GoResults]:
    """ Regenerate the previous results from JSON format. """
    if not previous:
        return None
    return Pfam2GoResults.from_json(previous, record)


def run_on_record(record: Record, results: Pfam2GoResults, options: ConfigType) -> Pfam2GoResults:
    """ Run the analysis, unless the previous results apply to the given record """
    if isinstance(results, Pfam2GoResults) and results.record_id == record.id:
        return results
    assert options.pfam2go
    return Pfam2GoResults(record.id, get_gos_for_pfams(record))
