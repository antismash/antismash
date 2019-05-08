# License: GNU Affero General Public License v3 or later
# A copy of GNU AGPL v3 should have been included in this software package in LICENSE.txt.

""" A collection of functions for running hmmsearch.
"""

from .base import execute, get_config, RunResult


def run_hmmpress(hmmfile: str) -> RunResult:
    """ Run hmmpress on a HMMer model, overwriting any previous generated files
        (e.g. '.h3i').

        Arguments:
            hmmfile: the path to the HMMer model

        Returns:
            a RunResult instance
    """
    return execute([get_config().executables.hmmpress, "-f", hmmfile])


def run_hmmpress_version() -> str:
    """ Get the version of the hmmpress """

    hmmpress = get_config().executables.hmmpress
    command = [
        hmmpress,
        "-h",
    ]

    help_text = execute(command).stdout
    if not help_text.startswith("# hmmpress"):
        msg = "unexpected output from hmmpress: %s, check path"
        raise RuntimeError(msg % hmmpress)

    version_line = help_text.split('\n')[1]
    return version_line.split()[2]
