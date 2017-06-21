#!/usr/bin/env python3
# vim: set fileencoding=utf-8 :
#
# License: GNU Affero General Public License v3 or later
# A copy of GNU AGPL v3 should have been included in this software package in LICENSE.txt.
"""Run the antiSMASH pipeline"""

import logging
import sys
import os

if sys.platform ==  ('win32') or sys.platform == ('darwin'):
    os.environ['EXEC'] = os.getcwd() + os.sep + "exec"
    os.environ['PYTHON'] = os.getcwd() + os.sep + "python"
    sys.path.append(os.sep.join([os.getcwd(), "python", "Lib", "site-packages"]))
    os.environ['PATH'] = os.pathsep + os.environ['PYTHON'] + os.pathsep + os.environ['PATH']
    os.environ['PATH'] = os.pathsep + os.environ['EXEC'] + os.pathsep + os.environ['PATH']


import antismash
from antismash.main import gather_modules
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

def main():
    parser = antismash.config.args.build_parser(from_config_file=True,
                                 modules=gather_modules(with_genefinding=True))

    #if --help, show help texts and exit
    if (list(set(["-h", "--help", "--help-showall"]) & set(sys.argv))):
        parser.print_help(None, "--help-showall" in sys.argv)
        return 0

    #Parse arguments, removing hyphens from the beginning of file names to avoid conflicts with argparse
    infile_extensions = ('.fasta', '.fas', '.fa', '.gb', '.gbk', '.emb', '.embl')
    sys.argv = [arg.replace("-","< > HYPHEN < >") if (arg.endswith(infile_extensions) and arg[0] == "-") else arg for arg in sys.argv]
    
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
        print(options.sequences)
        print("error: only one sequence file should be provided", file=sys.stderr)
        return 1
    if len(options.sequences) < 1:
        print("error: please specify a sequence file", file=sys.stderr)
        return 1
    sequence = options.sequences[0]
    del options.sequences
    config = antismash.config.args.Config(options)
   
    sequence = sequence.replace("< > HYPHEN < >","-")

    return antismash.run_antismash(sequence, config)


if __name__ == "__main__":
    sys.exit(main())
