# License: GNU Affero General Public License v3 or later
# A copy of GNU AGPL v3 should have been included in this software package in LICENSE.txt.

""" Lanthipeptides detection module

"""

import logging
from typing import Any, Dict, List, Optional

from antismash.common import path
from antismash.common.secmet import Record
from antismash.config import ConfigType
from antismash.config.args import ModuleArgs

from .config import get_config  # local config for fimo presence
from .specific_analysis import run_specific_analysis, LanthiResults
from .html_output import generate_html, will_handle

NAME = "lanthipeptides"
SHORT_DESCRIPTION = NAME.capitalize()


def get_arguments() -> ModuleArgs:
    """ Return the args for the lanthipeptide module """
    args = ModuleArgs('Advanced options', 'lanthi', enabled_by_default=True)
    return args


def check_options(_options: ConfigType) -> List:
    """ Lanthipeptide has no extra options, so there will be no conflicts """
    return []


def is_enabled(options: ConfigType) -> bool:
    """ Returns True if the lanthipeptide module is enabled """
    return options.lanthipeptides_enabled or not options.minimal


def regenerate_previous_results(results: Dict[str, Any], record: Record,
                                _options: ConfigType) -> Optional[LanthiResults]:
    """ Rebuilds the results from a prior run.

        Options aren't used here as the lanthipeptide module has no extra options.
    """
    if not results:
        return None
    regenned = LanthiResults.from_json(results, record)
    if not regenned:
        return None
    logging.debug("Reusing Lanthipeptide results: %d clusters contained %d total motifs",
                  len(regenned.clusters), sum(len(motifs) for motifs in regenned.motifs_by_locus.values()))
    return regenned


def check_prereqs() -> List[str]:
    """ Checks the prereqs for the lanthipeptide module.

        fimo is optional, having it available increases accuracy in the RODEO
        subsection
    """
    failure_messages = []
    for binary_name, optional in [('hmmpfam2', False), ('fimo', True)]:
        present = True
        if path.locate_executable(binary_name) is None:
            present = False
            if not optional:
                failure_messages.append("Failed to locate executable for %r" %
                                        binary_name)
        slot = '{}_present'.format(binary_name)
        conf = get_config()
        if hasattr(conf, slot):
            setattr(conf, slot, present)

    return failure_messages


def run_on_record(record: Record, results: LanthiResults, _options: ConfigType) -> LanthiResults:
    """ Runs the lanthipeptide analysis over the given record, if the existing
        results can't be reused.

        Options aren't used here as the lanthipeptide module has no extra options.
    """
    if isinstance(results, LanthiResults) and results.record_id == record.id:
        return results
    return run_specific_analysis(record)
