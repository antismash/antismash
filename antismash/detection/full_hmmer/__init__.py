# License: GNU Affero General Public License v3 or later
# A copy of GNU AGPL v3 should have been included in this software package in LICENSE.txt.

"""Full genome PFAM anotation"""

import logging
import os
from typing import List

from antismash.config import get_config
from antismash.common import fasta, path, subprocessing
from antismash.config.args import ModuleArgs

from .full_hmmer import generate_results, FullHmmerResults


NAME = "full_hmmer"
SHORT_DESCRIPTION = "Full genome PFAM anotation"


def get_arguments() -> ModuleArgs:
    """ Builds the module args """
    args = ModuleArgs('Full HMMer options', 'fullhmmer')
    args.add_analysis_toggle('fullhmmer',
                             dest='fullhmmer',
                             action='store_true',
                             default=False,
                             help="Run a whole-genome HMMer analysis.")
    return args


def check_options(options) -> List[str]:
    """ Never an issue as no args to check """
    return []


def is_enabled(options) -> bool:
    """  Uses the supplied options to determine if the module should be run """
    return options.fullhmmer


def check_prereqs() -> List[str]:
    "Check if all required applications are around"
    failure_messages = []
    for binary_name in ['hmmscan']:
        if not path.locate_executable(binary_name):
            failure_messages.append("Failed to locate file: %r" % binary_name)

    options = get_config()
    for file_name in ['Pfam-A.hmm', 'Pfam-A.hmm.h3f', 'Pfam-A.hmm.h3i',
                      'Pfam-A.hmm.h3m', 'Pfam-A.hmm.h3p']:
        if not path.locate_file(os.path.join(options.database_dir, 'pfam', file_name)):
            failure_messages.append("Failed to locate file: %r" % file_name)

    return failure_messages


def regenerate_previous_results(previous, record, _options) -> FullHmmerResults:
    """ Rebuild previous results """
    if not previous:
        return None
    return FullHmmerResults.from_json(previous, record)


def run_on_record(record, results, options) -> FullHmmerResults:
    "run hmmsearch against PFAM for all CDS features"
    if results:
        return results

    logging.info('Running whole-genome PFAM search')
    query_sequence = fasta.get_fasta_from_record(record)
    target_hmmfile = os.path.join(options.database_dir, 'pfam', 'Pfam-A.hmm')
    results = subprocessing.run_hmmscan(target_hmmfile, query_sequence)

    return generate_results(record, results, options)
