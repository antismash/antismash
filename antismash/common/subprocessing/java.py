# License: GNU Affero General Public License v3 or later
# A copy of GNU AGPL v3 should have been included in this software package in LICENSE.txt.

""" A collection of functions for running java.
"""

from .base import execute, get_config


def run_java_version() -> str:
    """ Get the version of the java binary """
    java = get_config().executables.java
    command = [
        java,
        "-version",
    ]

    version_string = execute(command).stderr
    if not version_string.startswith("openjdk") and not version_string.startswith("java"):
        msg = "unexpected output from java: %s, check path"
        raise RuntimeError(msg % java)
    # get rid of the non-version stuff in the output
    return version_string.split()[2].strip('"')
