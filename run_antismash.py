#!/usr/bin/env python3
# vim: set fileencoding=utf-8 :
#
# License: GNU Affero General Public License v3 or later
# A copy of GNU AGPL v3 should have been included in this software package in LICENSE.txt.
"""Run the antiSMASH pipeline"""

import logging
import sys
import os

import antismash
from antismash.common.subprocessing import execute


def get_git_version():
    """Get the sha1 of the current git version"""
    args = ['git', 'rev-parse', '--short', 'HEAD']
    try:
        return execute(args).stdout.strip()
    except OSError:
        pass

    return ""

def get_version():
    """Get the current version string"""
    version = antismash.__version__
    git_version = get_git_version()
    if git_version != '':
        version += "-%s" % git_version

    return version

def main(args):
    all_modules = antismash.get_detection_modules() + antismash.get_analysis_modules()
    parser = antismash.config.args.build_parser(from_config_file=True, modules=all_modules)

    #if --help, show help texts and exit
    if (list(set(["-h", "--help", "--help-showall"]) & set(args))):
        parser.print_help(None, "--help-showall" in args)
        return 0

    #Parse arguments, removing hyphens from the beginning of file names to avoid conflicts with argparse
    infile_extensions = ('.fasta', '.fas', '.fa', '.gb', '.gbk', '.emb', '.embl')
    args = [arg.replace("-","< > HYPHEN < >") if (arg.endswith(infile_extensions) and arg[0] == "-") else arg for arg in args]

    try:
        options = parser.parse_args(["@config_test"] + sys.argv[1:])
    except SystemExit:
        # note: logging isn't set up as usual here because it relies on config
        logging.error("option generation exited early")
        raise

    #if -V, show version text and exit
    if options.version:
        print("antiSMASH %s" % get_version())
        return 0

    if len(options.sequences) > 1:
        parser.error("Only one sequence file should be provided")
        return 1
    if len(options.sequences) < 1 and not options.reuse_results:
        parser.error("One of an input file or --reuse-results must be specified")
        return 1
    if options.sequences:
        sequence = options.sequences[0]
        del options.sequences
    else:
        sequence = ""

    # if not supplied, set the output directory to be the sequence name
    # can't be done in argparse because parsing interacting args is a bad idea
    if not options.output_dir:
        options.output_dir = os.path.abspath(os.path.splitext(os.path.basename(sequence))[0])

    config = antismash.config.args.Config(options)

    sequence = sequence.replace("< > HYPHEN < >","-")

    return antismash.run_antismash(sequence, config)


if __name__ == "__main__":
    sys.exit(main(sys.argv))
