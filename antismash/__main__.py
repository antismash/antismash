# vim: set fileencoding=utf-8 :
#
# License: GNU Affero General Public License v3 or later
# A copy of GNU AGPL v3 should have been included in this software package in LICENSE.txt.
"""Run the antiSMASH pipeline"""

import os
import sys
from typing import List, Optional
# wrap the antismash imports in a try/except as the import is lengthy and SIGINT
# is noisy when it happens during the import process
try:

    import antismash
    from antismash.common.path import (
        changed_directory,
        get_full_path,
        locate_file,
    )
    from antismash.common.subprocessing import execute
except KeyboardInterrupt:
    sys.exit(2)

GIT_VERSION_FALLBACK_FILENAME = get_full_path(__file__, "git_hash")


def get_git_version(fallback_filename: Optional[str] = GIT_VERSION_FALLBACK_FILENAME) -> str:
    """Get the sha1 of the current git version"""
    git_version = ""
    try:
        with changed_directory(os.path.dirname(__file__)):
            version_cmd = execute(['git', 'rev-parse', '--short', 'HEAD'])
            status_cmd = execute(['git', 'status', '--porcelain'])
        if version_cmd.successful() and status_cmd.successful():
            git_version = version_cmd.stdout.strip()
            changes = status_cmd.stdout.splitlines()
            if changes:
                git_version += "(changed)"
    except OSError:
        pass
    if git_version == "" and fallback_filename:
        if locate_file(fallback_filename, silent=True):
            with open(fallback_filename, "rt", encoding="utf-8") as handle:
                git_version = handle.read().strip()
    return git_version


def get_version() -> str:
    """Get the current version string"""
    # if already set, use that version rather than redo the git work
    if antismash.config.args.ANTISMASH_VERSION:
        return antismash.config.args.ANTISMASH_VERSION
    version = antismash.__version__
    git_version = get_git_version()
    if git_version:
        version += f"-{git_version}"
    return version


# set arg version to avoid cyclic imports
antismash.config.args.ANTISMASH_VERSION = get_version()


def main(args: List[str], *, branding_override: str = "", version_override: str = "") -> int:
    """ The entrypoint of antiSMASH as if it was on the command line

        Arguments:
            args: a list of args as would be given on the command line
                    e.g. ["inputfile", "--minimal", "--enable-nrps_pks"]
            branding_override: a name to use instead of antiSMASH, when used as a library
            branding_override: a version to use instead of antiSMASH's version, when used as a library

        Returns:
            zero if successful, non-zero otherwise

    """
    all_modules = antismash.get_all_modules()
    parser = antismash.config.args.build_parser(from_config_file=True, modules=all_modules)

    if branding_override or version_override:
        if not branding_override and version_override:
            raise ValueError("both branding and version must be overridden, if either is supplied")
        assert branding_override != "antiSMASH"
        antismash.config.args.ANTISMASH_VERSION = (
            f"{ version_override }"
            f" (based on antiSMASH { get_version() })"
        )

    # if --help, show help texts and exit
    if set(args).intersection({"-h", "--help", "--help-showall"}):
        parser.print_help(show_all="--help-showall" in args)
        return 0

    options = antismash.config.build_config(args, parser=parser)
    # while the branding can be set in the default configs, ensure it's set here too
    if branding_override and options.branding != branding_override:
        antismash.config.update_config({"branding": branding_override})

    assert options.branding, repr(options.branding)

    if options.write_config_file:
        parser.write_to_config_file(options.write_config_file, options.__dict__)
        return 0

    # if -V, show version text and exit
    if options.version:
        print(f"{options.branding} {get_version()}")
        return 0

    if len(options.sequences) > 1:
        parser.error("Only one sequence file should be provided")
    if len(options.sequences) < 1 and not options.reuse_results \
            and not options.check_prereqs_only and not options.list_plugins:
        parser.error("One of an input file or --reuse-results must be specified")
    if options.sequences and options.reuse_results:
        parser.error("Provide a sequence file or results to reuse, not both.")
    if options.sequences:
        sequence = options.sequences[0]
        if not os.path.exists(sequence):
            parser.error(f"Input file does not exist: {sequence}")
        if not os.path.isfile(sequence):
            parser.error(f"input {sequence} is not a file")
    else:
        sequence = ""

    if options.reuse_results and not os.path.exists(options.reuse_results):
        parser.error(f"Input file does not exist: {options.reuse_results}")

    if os.sep in options.output_basename:
        parser.error(f"Output basename cannot contain a path separator character {os.sep}")

    options.version = get_version()

    try:
        return antismash.run_antismash(sequence, options)
    except antismash.common.errors.AntismashInputError as err:
        if not str(err):
            raise
        print("ERROR:", str(err), file=sys.stderr)
        return 1


def entrypoint() -> None:
    """This is needed for the script generated by setuptools."""
    try:
        sys.exit(main(sys.argv[1:]))
    except KeyboardInterrupt:
        print("ERROR: Interrupted", file=sys.stderr)
        sys.exit(2)


if __name__ == '__main__':
    entrypoint()
