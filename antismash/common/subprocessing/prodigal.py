# License: GNU Affero General Public License v3 or later
# A copy of GNU AGPL v3 should have been included in this software package in LICENSE.txt.

""" A collection of functions for running prodigal.
"""

from .base import execute, get_config


def run_prodigal_version() -> str:
    """ Get the version of the prodigal binary """
    prodigal = get_config().executables.prodigal
    command = [
        prodigal,
        "-v",
    ]

    version_string = execute(command).stderr.strip()
    if not version_string.startswith("Prodigal"):
        msg = "unexpected output from prodigal: %s, check path"
        raise RuntimeError(msg % prodigal)
    # get rid of the non-version stuff in the output
    return version_string.split()[1][:-1]
