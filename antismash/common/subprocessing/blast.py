# License: GNU Affero General Public License v3 or later
# A copy of GNU AGPL v3 should have been included in this software package in LICENSE.txt.

""" A collection of functions for running ncbi-blast+ commands.
"""

from tempfile import NamedTemporaryFile
from typing import Any, IO, List

from .base import execute, get_config, RunResult, SearchIO

NCBI_OVERRIDES = {
    # NCBI tools have added usage reporting, which can add minutes to simple
    # things like fetching the version, this disables that reporting
    "BLAST_USAGE_REPORT": "false",
}


def run_blastp(target_blastp_database: str, query_sequence: str,
               opts: List[str] = None, results_file: str = None
               ) -> List[SearchIO._model.query.QueryResult]:
    """ Runs blastp over a single sequence against a database and returns the
        results as parsed by Bio.SearchIO.

        Arguments:
            target_blastp_database: the blastp database to compare to
            query_sequence: the sequence being compared
            opts: a list of extra arguments to pass to blastp, or None
            results_file: a path to keep a copy of blastp results in, if provided

        Returns:
            a list of QueryResults as parsed from blast output by SearchIO
    """
    if not query_sequence:
        raise ValueError("Cannot run blastp on empty sequence")

    config = get_config()
    command = [
        config.executables.blastp,
        "-num_threads", str(config.cpus),
        "-db", target_blastp_database,
        "-outfmt", "6",  # use tabular format for biopython's parsing
    ]

    if opts is not None:
        command.extend(opts)

    if results_file is not None:
        handle: IO[Any] = open(results_file)  # pylint: disable=consider-using-with
    else:
        handle = NamedTemporaryFile()  # pylint: disable=consider-using-with

    result = execute(command, stdin=query_sequence, stdout=handle, environment_overrides=NCBI_OVERRIDES)
    if not result.successful():
        raise RuntimeError('blastp returned %d: %r while scanning %r' % (
                           result.return_code, result.stderr.replace("\n", ""),
                           query_sequence[:100]))
    filename = results_file or handle.name
    with open(filename) as read_handle:
        parsed = list(SearchIO.parse(read_handle, 'blast-tab'))
    handle.close()
    return parsed


def run_blastp_version() -> str:
    """ Get the version of the blastp binary """
    blastp = get_config().executables.blastp
    command = [
        blastp,
        "-version",
    ]

    version_string = execute(command, environment_overrides=NCBI_OVERRIDES).stdout
    if not version_string.startswith("blastp: "):
        msg = "unexpected output from blastp: %s, check path"
        raise RuntimeError(msg % blastp)
    # get rid of the non-version stuff in the output
    return version_string.split('\n')[0][8:]


def run_makeblastdb(filename: str, target: str = "", dbtype: str = "prot") -> RunResult:
    """ Runs makeblastdb on the inputs to create a blast database, creating
        additional files with the given target name.

        Arguments:
            filename: the input filename
            target: the output name (including directory) for the created database
                    (defaults to input name)
            dbtype: the type of database (defaults to 'prot')

        Returns:
            a subprocessing.RunResult instance
    """
    if not target:
        target = filename
    command = [
        get_config().executables.makeblastdb,
        "-in", filename,
        "-out", target,
        "-dbtype", dbtype,
    ]
    result = execute(command, environment_overrides=NCBI_OVERRIDES)
    if not result.successful():
        raise RuntimeError(f"makeblastdb failed to run: {command} -> {result.stderr.splitlines()[-1]}")
    return result


def run_makeblastdb_version() -> str:
    """ Get the version of the makeblastdb binary """
    makeblastdb = get_config().executables.makeblastdb
    command = [
        makeblastdb,
        "-version",
    ]

    version_string = execute(command, environment_overrides=NCBI_OVERRIDES).stdout
    if not version_string.startswith("makeblastdb: "):
        msg = "unexpected output from makeblastdb: %s, check path"
        raise RuntimeError(msg % makeblastdb)
    # get rid of the non-version stuff in the output
    return version_string.split('\n')[0][13:]
