# License: GNU Affero General Public License v3 or later
# A copy of GNU AGPL v3 should have been included in this software package in LICENSE.txt.

""" General HMM detection module for detecting specific domains and defining
    clusters based on domains detected
"""

import os
from typing import List

import antismash.common.path as path
from antismash.common.subprocessing import run_hmmpress

from antismash.config.args import ModuleArgs

from antismash.detection.hmm_detection import rule_parser
from antismash.detection.hmm_detection.hmm_detection import detect_signature_genes
from antismash.detection.hmm_detection.rule_parser import Parser
from antismash.detection.hmm_detection.signatures import get_signature_profiles

NAME = "hmmdetection"
SHORT_DESCRIPTION = "Cluster detection with HMM"


def get_supported_cluster_types() -> List[str]:
    """ Returns a list of all cluster types for which there are rules
    """
    with open(path.get_full_path(__file__, 'cluster_rules.txt'), "r") as rulefile:
        rules = rule_parser.Parser("".join(rulefile.readlines())).rules
        clustertypes = [rule.name for rule in rules]
    return clustertypes


def get_arguments() -> ModuleArgs:
    """ Constructs commandline arguments and options for this module
    """
    args = ModuleArgs('Advanced options', '', override_safeties=True)
    cluster_types = get_supported_cluster_types()
    args.add_option('--enable',
                    metavar="TYPES",
                    dest='enabled_cluster_types',
                    type=lambda x: x.split(","),
                    default=cluster_types,
                    help=("Select sec. met. cluster types to search for. "
                          " E.g. --enable t1pks,nrps,other"))
    return args


def check_options(options) -> List[str]:
    """ Checks the options to see if there are any issues before
        running any analyses
    """
    available = set(get_supported_cluster_types())
    enabled = options.enabled_cluster_types
    if isinstance(enabled, list):
        requested = set(enabled)
    else:
        requested = set(enabled.replace(";", ",").split(","))

    unavailable = requested - available

    errors = []
    if unavailable:
        errors.append("The following cluster types are unavailable:")
        errors.extend(" "+cluster for cluster in unavailable)
        errors.append("Available types are:")
        errors.extend(" "+cluster for cluster in sorted(available))
    return errors


def is_enabled(_options) -> bool:
    """  Uses the supplied options to determine if the module should be run
    """
    # in this case, yes, always
    return True


def regenerate_previous_results(results, record, options) -> None:  # pylint: disable=unused-argument
    """ Regenerate previous results. """
    # always rerun hmmdetection  # TODO: should clusters be kept?
    return None


def run_on_record(record, options) -> None:
    """ Runs hmm_detection on the provided record.
    """
    return detect_signature_genes(record, options)


def check_prereqs() -> List[str]:
    """ Check that prereqs are satisfied. hmmpress is only required if the
        databases have not yet been generated.
    """
    failure_messages = []
    for binary_name, optional in [('hmmsearch', False), ('hmmpress', False)]:
        if path.locate_executable(binary_name) is None and not optional:
            failure_messages.append("Failed to locate executable for %r" %
                                    binary_name)

    profiles = None
    # Check that hmmdetails.txt is readable and well-formatted
    try:
        profiles = get_signature_profiles()
    except ValueError as err:
        failure_messages.append(str(err))

    # the path to the markov model
    hmm = path.get_full_path(__file__, 'data', 'bgc_seeds.hmm')
    hmm_files = [os.path.join("data", sig.hmm_file) for sig in profiles]
    if path.locate_file(hmm) is None:
        # try to generate file from all specified profiles in hmmdetails
        try:
            with open(hmm, 'w') as all_hmms_handle:
                for hmm_file in hmm_files:
                    with open(path.get_full_path(__file__, hmm_file), 'r') as handle:
                        all_hmms_handle.write(handle.read())
        except OSError:
            failure_messages.append('Failed to generate file {!r}'.format(hmm))

    # if previous steps have failed, the remainder will too, so don't try
    if failure_messages:
        return failure_messages

    # Check that cluster_rules.txt is readable and well-formatted
    try:
        with open(path.get_full_path(__file__, "cluster_rules.txt")) as rules:
            parser = Parser("".join(rules.readlines()))
        if not parser.rules:
            failure_messages.append("No rules contained in cluster_rules.txt")
    except ValueError as err:
        failure_messages.append(str(err))

    binary_extensions = ['.h3f', '.h3i', '.h3m', '.h3p']
    for ext in binary_extensions:
        binary = "{}{}".format(hmm, ext)
        if path.locate_file(binary) is None:
            result = run_hmmpress(hmm)
            if not result.succesful():
                failure_messages.append('Failed to hmmpress {!r}: {}'.format(hmm, result.stderr))
            break

    return failure_messages
