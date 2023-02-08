# License: GNU Affero General Public License v3 or later
# A copy of GNU AGPL v3 should have been included in this software package in LICENSE.txt.

""" Detection of NRPS/PKS domains within genes
"""

import logging
from typing import Any, Dict, List, Optional

from antismash.common import hmmer, path
from antismash.common.secmet import Record
from antismash.config import ConfigType, get_config
from antismash.config.args import ModuleArgs
from antismash.detection import DetectionStage

from .domain_drawing import generate_html, generate_js_domains, will_handle
from .domain_identification import generate_domains, get_database_path, NRPSPKSDomains
from .modular_domain import ModularDomain

NAME = "nrps_pks_domains"
SHORT_DESCRIPTION = "NRPS/PKS domain identification"
DETECTION_STAGE = DetectionStage.PER_AREA


def get_arguments() -> ModuleArgs:
    """ Constructs commandline arguments and options for this module
    """
    args = ModuleArgs('Advanced options', '', override_safeties=True)
    return args


def check_options(_options: ConfigType) -> List[str]:
    """ Checks the options to see if there are any issues before
        running any analyses
    """
    return []


def is_enabled(_options: ConfigType) -> bool:
    """  Uses the supplied options to determine if the module should be run
    """
    # in this case, yes, always
    return True


def regenerate_previous_results(results: Dict[str, Any], record: Record,
                                _options: ConfigType) -> Optional[NRPSPKSDomains]:
    """ Reconstruct NRPS/PKS domain detection results from a JSON format """
    return NRPSPKSDomains.from_json(results, record)


def run_on_record(record: Record, previous_results: Optional[NRPSPKSDomains],
                  _options: ConfigType) -> NRPSPKSDomains:
    """ Find and mark NRPS and PKS domains within CDSFeatures contained in clusters """
    logging.debug('Marking NRPS/PKS CDS features and domains in clusters')
    if previous_results:
        assert isinstance(previous_results, NRPSPKSDomains)
        return previous_results
    results = generate_domains(record)
    return results


def prepare_data(logging_only: bool = False) -> List[str]:
    """ Ensures packaged data is fully prepared

        Arguments:
            logging_only: whether to return error messages instead of raising exceptions

        Returns:
            a list of error messages (only if logging_only is True)
    """
    failure_messages = []
    for model in ['abmotifs.hmm', 'dockingdomains.hmm', 'ksdomains.hmm', 'nrpspksdomains.hmm']:
        full_path = path.get_full_path(__file__, "data", model)
        failure_messages.extend(hmmer.ensure_database_pressed(full_path, return_not_raise=logging_only))
    options = get_config()
    if "mounted_at_runtime" in options.database_dir:
        return failure_messages
    for subdir, filename in [("transATor", "transATor.hmm")]:
        try:
            full_path = get_database_path(subdir, filename)
        except ValueError as err:
            if logging_only:
                failure_messages.append(str(err))
                continue
            raise
        failure_messages.extend(hmmer.ensure_database_pressed(full_path, return_not_raise=logging_only))
    return failure_messages


def check_prereqs(options: ConfigType) -> List[str]:
    """ Ensures the various required hmmer profile files exist """
    failure_messages = []
    for binary_name in ['hmmscan', 'hmmpress']:
        if binary_name not in options.executables:
            failure_messages.append(f"Failed to locate executable for {binary_name!r}")

    # skip trying to press files if there's missing binaries
    if failure_messages:
        return failure_messages

    failure_messages.extend(prepare_data(logging_only=True))

    return failure_messages
