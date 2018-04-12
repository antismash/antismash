#!/usr/bin/env python3
# vim: set fileencoding=utf-8 :
#
# License: GNU Affero General Public License v3 or later
# A copy of GNU AGPL v3 should have been included in this software package in LICENSE.txt.
"""Run the antiSMASH pipeline"""

import sys
from typing import List

import antismash
from antismash.common.subprocessing import execute
from antismash.config import ConfigType


def get_git_version() -> str:
    """Get the sha1 of the current git version"""
    args = ['git', 'rev-parse', '--short', 'HEAD']
    try:
        return execute(args).stdout.strip()
    except OSError:
        pass
    return ""


def get_version() -> str:
    """Get the current version string"""
    version = antismash.__version__
    git_version = get_git_version()
    if git_version:
        version += "-%s" % git_version

    return version


def main(args: List[str]) -> int:
    all_modules = antismash.get_detection_modules() + antismash.get_analysis_modules()
    parser = antismash.config.args.build_parser(from_config_file=True, modules=all_modules)

    # if --help, show help texts and exit
    if set(args).intersection({"-h", "--help", "--help-showall"}):
        parser.print_help(None, "--help-showall" in args)
        return 0

    options = antismash.config.build_config(args, parser=parser)

    if options.write_config_file:
        parser.write_to_config_file(options.write_config_file, options.__dict__)
        return 0

    # if -V, show version text and exit
    if options.version:
        print("antiSMASH %s" % get_version())
        return 0

    if len(options.sequences) > 1:
        parser.error("Only one sequence file should be provided")
        return 1
    if len(options.sequences) < 1 and not options.reuse_results \
            and not options.check_prereqs_only and not options.list_plugins:
        parser.error("One of an input file or --reuse-results must be specified")
        return 1
    if options.sequences and options.reuse_results:
        parser.error("Provide a sequence file or results to reuse, not both.")
        return 1
    if options.sequences:
        sequence = options.sequences[0]
        options.__dict__.pop("sequences")
    else:
        sequence = ""

    options.version = get_version()

    return antismash.run_antismash(sequence, options)


if __name__ == "__main__":
    sys.exit(main(sys.argv[1:]))
