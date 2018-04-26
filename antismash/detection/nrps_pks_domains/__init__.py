# License: GNU Affero General Public License v3 or later
# A copy of GNU AGPL v3 should have been included in this software package in LICENSE.txt.

""" Detection of NRPS/PKS domains within genes
"""

import logging
from typing import Any, Dict, List, Optional

from antismash.common import path, subprocessing
from antismash.common.secmet import Record
from antismash.config import ConfigType
from antismash.config.args import ModuleArgs

from .domain_identification import generate_domains, NRPSPKSDomains

NAME = "nrps_pks_domains"
SHORT_DESCRIPTION = "NRPS/PKS domain identification"


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


def regenerate_previous_results(results: Dict[str, Any], record: Record, _options: ConfigType) -> NRPSPKSDomains:
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
    results.add_to_record(record)
    return results


def check_prereqs() -> List[str]:
    """ Ensures the various required hmmer profile files exist """
    failure_messages = []
    for binary_name, optional in [('hmmscan', False), ('hmmpress', False)]:
        if path.locate_executable(binary_name) is None and not optional:
            failure_messages.append("Failed to locate executable for %r" %
                                    binary_name)

    markov_models = [path.get_full_path(__file__, 'data', filename) for filename in [
                                'abmotifs.hmm', 'dockingdomains.hmm',
                                'ksdomains.hmm', 'nrpspksdomains.hmm']]

    binary_extensions = ['.h3f', '.h3i', '.h3m', '.h3p']

    for hmm in markov_models:
        if path.locate_file(hmm) is None:
            failure_messages.append("Failed to locate file %r" % hmm)
            continue
        for ext in binary_extensions:
            binary = "{}{}".format(hmm, ext)
            if path.locate_file(binary) is None:
                result = subprocessing.run_hmmpress(hmm)
                if not result.successful():
                    failure_messages.append('Failed to hmmpress {!r}: {}'.format(hmm, result.stderr))
                break

    return failure_messages
