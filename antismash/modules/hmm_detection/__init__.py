# License: GNU Affero General Public License v3 or later
# A copy of GNU AGPL v3 should have been included in this software package in LICENSE.txt.

"""General HMM detection module

"""
import argparse
import logging
import os

import antismash.common.deprecated as utils
import antismash.common.path as path
from antismash import config
from antismash.common.deprecated import SeqFeature, FeatureLocation
from antismash.common.subprocessing import run_hmmsearch, run_hmmpress

from antismash.config.args import ModuleArgs

from antismash.modules.hmm_detection import rule_parser
from antismash.modules.hmm_detection.hmm_detection import detect_signature_genes
from antismash.modules.hmm_detection.rule_parser import Parser
from antismash.modules.hmm_detection.signatures import get_signature_profiles

NAME = "hmmdetection" # used when building commandline args for dis-/enabling
SHORT_DESCRIPTION = "Cluster detection with HMM"

def get_supported_cluster_types():
    """ Returns a list of all cluster types for which there are rules
    """
    with open(path.get_full_path(__file__, 'cluster_rules.txt'), "r") as rulefile:
        rules = rule_parser.Parser(rulefile).rules
        clustertypes = [rule.name for rule in rules]
    return clustertypes

def get_arguments():
    """ Constructs commandline arguments and options for this module
    """
    args = ModuleArgs('Advanced options', '', override_safeties=True)
    cluster_types = get_supported_cluster_types()
    args.add_argument('--enable',
                       metavar="TYPES",
                       dest='enabled_cluster_types',
                       type=str,
                       default=cluster_types,
                       help="Select sec. met. cluster types to search for. E.g. --enable t1pks,nrps,other")
    return args

def check_options(options):
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

def is_enabled(options):
    """  Uses the supplied options to determine if the module should be run
    """
    # in this case, yes, always
    return True

def run_on_record(seq_record, options):
    return detect_signature_genes(seq_record, options)


def check_prereqs():
    failure_messages = []
    for binary_name, optional in [('hmmsearch', False), ('hmmpress', False)]:
        if path.locate_executable(binary_name) is None and not optional:
            failure_messages.append("Failed to locate executable for %r" %
                                    binary_name)

    profiles = None
    # Check that hmmdetails.txt is readable and well-formatted
    try:
        profiles = get_signature_profiles()
    except ValueError as e:
        failure_messages.append(str(e))


    # the path to the markov model
    hmm = path.get_full_path(__file__, 'data/bgc_seeds.hmm')
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

    # if previous steps have failed, the remaineder will too, so don't try
    if failure_messages:
        return failure_messages

    # Check that cluster_rules.txt is readable and well-formatted
    try:
        with open(path.get_full_path(__file__, "cluster_rules.txt")) as rules:
            parser = Parser(rules.readlines())
        if not parser.rules:
            failure_messages.append("No rules contained in cluster_rules.txt")
    except ValueError as e:
        failure_messages.append(str(e))


    binary_extensions = ['.h3f', '.h3i', '.h3m', '.h3p']
    for ext in binary_extensions:
        binary = "{}{}".format(hmm, ext)
        if path.locate_file(binary) is None:
            _, err, retcode = run_hmmpress(hmm)
            if retcode != 0:
                failure_messages.append('Failed to hmmpress {!r}: {!r}'.format(hmm, err))
            break

    return failure_messages

__all__ = ["check_prereqs", "detect_signature_genes"]
