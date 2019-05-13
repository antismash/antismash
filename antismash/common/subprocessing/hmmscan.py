# License: GNU Affero General Public License v3 or later
# A copy of GNU AGPL v3 should have been included in this software package in LICENSE.txt.

""" A collection of functions for running hmmscan.
"""

from io import StringIO
from typing import List

from .base import execute, get_config, SearchIO


def run_hmmscan(target_hmmfile: str, query_sequence: str, opts: List[str] = None,
                results_file: str = None) -> List[SearchIO._model.query.QueryResult]:  # pylint: disable=protected-access
    """ Runs hmmscan on the inputs and return a list of QueryResults

        Arguments:
            target_hmmfile: the path to a HMM file to use in scanning
            query_sequence: a string containing input sequences in fasta format
            opts: a list of extra arguments to pass to hmmscan, or None
            results_file: a path to keep a copy of hmmscan results in, if provided

        Returns:
            a list of QueryResults as parsed from hmmscan output by SearchIO

    """
    if not query_sequence:
        raise ValueError("Cannot run hmmscan on empty sequence")

    config = get_config()
    command = [config.executables.hmmscan, "--cpu", str(config.cpus), "--nobias"]

    # Only run multithreaded when the binary supports it
    if " --cpu " not in run_hmmscan_help():
        command = command[0:1] + command[3:]

    if opts is not None:
        command.extend(opts)
    command.extend([target_hmmfile, '-'])
    result = execute(command, stdin=query_sequence)
    if not result.successful():
        raise RuntimeError('hmmscan returned %d: %r while scanning %r' % (
                           result.return_code, result.stderr[-100:].replace("\n", ""),
                           query_sequence[:100]))
    if results_file is not None:
        with open(results_file, 'w') as handle:
            handle.write(result.stdout)

    return list(SearchIO.parse(StringIO(result.stdout), 'hmmer3-text'))


def run_hmmscan_help() -> str:
    """ Get the help output of hmmscan """
    # cache results
    help_text = getattr(run_hmmscan_help, 'help_text', '')
    if help_text:
        return help_text

    hmmscan = get_config().executables.hmmscan
    command = [
        hmmscan,
        "-h",
    ]

    help_text = execute(command).stdout
    if not help_text.startswith("# hmmscan"):
        msg = "unexpected output from hmmscan: %s, check path"
        raise RuntimeError(msg % hmmscan)

    setattr(run_hmmscan_help, 'help_text', help_text)
    return help_text


def run_hmmscan_version() -> str:
    """ Get the version of the hmmscan """
    version_line = run_hmmscan_help().split('\n')[1]
    return version_line.split()[2]
