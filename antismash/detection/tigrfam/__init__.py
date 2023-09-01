# License: GNU Affero General Public License v3 or later
# A copy of GNU AGPL v3 should have been included in this software package in LICENSE.txt.

""" TIGRFam anotation for only clusters """

import logging
import os
from typing import Any, Dict, List, Optional

from antismash.common import hmmer, path
from antismash.common.secmet import Record
from antismash.config import ConfigType
from antismash.config.args import ModuleArgs
from antismash.detection import DetectionStage

from .tigr_results import TIGRFamResults

NAME = "tigrfam"
SHORT_DESCRIPTION = "Cluster TIGRFam anotation"
DETECTION_STAGE = DetectionStage.PER_AREA

MIN_SCORE = 0.
MAX_EVALUE = 0.01


def get_arguments() -> ModuleArgs:
    """ Builds the module args """
    args = ModuleArgs('TIGRFam options', 'tigrfam')
    args.add_analysis_toggle('tigrfam',
                             dest='tigrfam',
                             action='store_true',
                             default=False,
                             help="Annotate clusters using TIGRFam profiles.")
    return args


def is_enabled(options: ConfigType) -> bool:
    """  Uses the supplied options to determine if the module should be run """
    return options.tigrfam


def check_prereqs(options: ConfigType) -> List[str]:
    """ Ensure at least one database exists and is valid """
    failure_messages = []
    for binary_name in ['hmmscan']:
        if binary_name not in options.executables:
            failure_messages.append(f"Failed to locate executable: {binary_name!r}")

    # account for database directories mounted into docker containers
    if "mounted_at_runtime" in options.database_dir:
        return failure_messages

    tigr_db = os.path.join(options.database_dir, "tigrfam", "TIGRFam.hmm")
    if not path.locate_file(tigr_db):
        failure_messages.append(f"Failed to locate TIGRFam db in {os.path.join(options.database_dir, 'tigrfam')}")

    failure_messages.extend(hmmer.ensure_database_pressed(tigr_db, return_not_raise=True))

    return failure_messages


def check_options(_options: ConfigType) -> List[str]:
    """ Check the requested PFAM database exists """
    # No options to check
    return []


def regenerate_previous_results(previous: Dict[str, Any], record: Record,
                                _options: ConfigType) -> Optional[TIGRFamResults]:
    """ Rebuild previous results """
    if not previous:
        return None
    results = TIGRFamResults.from_json(previous, record)
    if not results:
        return None
    if results.score > MIN_SCORE or results.evalue < MAX_EVALUE:
        # new values too lenient, discard resuts
        return None
    return results.refilter(MAX_EVALUE, MIN_SCORE)


def run_on_record(record: Record, results: Optional[TIGRFamResults],
                  options: ConfigType) -> TIGRFamResults:
    """ Run hmmsearch against TIGRFam for all CDS features within the record """

    logging.info('Running TIGRFam search')

    if results:
        return results

    features = record.get_cds_features_within_regions()
    tigr_db = os.path.join(options.database_dir, "tigrfam", "TIGRFam.hmm")
    hmmer_results = hmmer.run_hmmer(record, features, MAX_EVALUE, MIN_SCORE, tigr_db,
                                    "tigrfam", filter_overlapping=False)
    return TIGRFamResults.from_hmmer_results(hmmer_results)
