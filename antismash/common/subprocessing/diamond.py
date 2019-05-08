# License: GNU Affero General Public License v3 or later
# A copy of GNU AGPL v3 should have been included in this software package in LICENSE.txt.

""" A collection of functions for running diamond.
"""


from typing import List, Optional

from helperlibs.wrappers.io import TemporaryDirectory

from .base import execute, get_config, RunResult


def run_diamond(subcommand: str,
                opts: Optional[List[str]] = None) -> RunResult:
    """ Run a diamond subcommand, possibly with further options.

        Arguments:
            subcommand: the diamond subcommand to run
            opts: a list of additional argument strings to pass to diamond

        Returns:
            RunResult of running diamond
    """
    config = get_config()
    if not config.executables.diamond:
        raise RuntimeError("no available diamond executable")
    with TemporaryDirectory() as temp_dir:
        params = [
            config.executables.diamond,
            subcommand,
            "--threads", str(config.cpus),
            "--tmpdir", temp_dir,
        ]

        if opts:
            params.extend(opts)

        result = execute(params)
        if not result.successful():
            raise RuntimeError("diamond failed to run: %s -> %s" % (subcommand, result.stderr[-100:]))
    return result


def run_diamond_search(query_file: str, database_file: str, mode: str = "blastp",
                       opts: Optional[List[str]] = None) -> str:
    """ Runs diamond, comparing the given query to the given database

        Arguments:
            query_file: the path of query sequence file
            database_file: the path of the database to compare to
            mode: the mode to use (defaults to blastp)
            opts: any extra options to pass to diamond

        Returns:
            the output from running diamond
    """

    args = [
        "--db", database_file,
        "--query", query_file,
    ]

    if opts:
        args.extend(opts)

    return run_diamond(mode, args).stdout


def run_diamond_makedb(database_file: str, sequence_file: str) -> RunResult:
    """ Generate a new diamond database from a fasta file.Optional

        Arguments:
            database_file: the path to the database to generate
            sequence_file: the path to a protein FASTA file to generate the database from

        Returns:
            the RunResult running diamond
    """
    args = [
        "--db", database_file,
        "--in", sequence_file,
    ]

    return run_diamond("makedb", args)


def run_diamond_version() -> str:
    """ Get the version of the diamond executable

        Returns:
            The numeric part of "diamond version"
    """

    version_string = run_diamond("version").stdout
    if not version_string.startswith("diamond version "):
        msg = "unexpected output from diamond-executable: %s, check path"
        raise RuntimeError(msg % get_config().executables.diamond)
    # Get rid of the "diamond version" prefix
    return version_string[16:].strip()
