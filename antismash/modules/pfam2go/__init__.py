# License: GNU Affero General Public License v3 or later
# A copy of GNU AGPL v3 should have been included in this software package in LICENSE.txt.

"""Map Pfam domains to Gene Ontology terms"""

from typing import Any, Dict, List, Optional

from antismash.config.args import ModuleArgs
from antismash.modules.pfam2go.pfam2go import parse_all_mappings, build_at_the_end, build_as_i_go,\
    construct_gene_ontologies, GeneOntology, GeneOntologies, Pfam2GoResults, get_gos_for_pfams

NAME = "pfam2go"
SHORT_DESCRIPTION = "Pfam domain to GO mapping"


def get_arguments() -> ModuleArgs:
    """ Build and return arguments. No extra options beyond a switch to enable """
    args = ModuleArgs('Additional analysis', 'pfam2go')
    args.add_analysis_toggle('--pfam2go',
                             dest='pfam2go',
                             action='store_true',
                             default=False,
                             help="Run Pfam to Gene Ontology mapping module.")
    return args


def check_options(_options) -> List[str]:
    """ Checks options for conflicts.
        No extra options, so they can't have conflicts.
    """
    return []


def check_prereqs() -> List[str]:
    """Check for prerequisites"""
    # Needs Pfam IDs, how to check?
    return []


def is_enabled(options) -> bool:
    """ Should the module be run with these options """
    return options.pfam2go


def regenerate_previous_results(previous: Dict[str, Any], record, _options) -> Optional[Pfam2GoResults]:
    """ Regenerate the previous results from JSON format. """
    if not previous:
        return None
    return Pfam2GoResults.from_json(previous, record)


def run_on_record(record, results: Pfam2GoResults, options) -> Pfam2GoResults:
    """ Run the analysis, unless the previous results apply to the given record """
    if isinstance(results, Pfam2GoResults) and results.record_id == record.id:
        return results
    #  otherwise, extract Pfam IDs, do pfam to GO mapping with these?
    #return detect(record, options)
    return Pfam2GoResults(record.id, get_gos_for_pfams(record))