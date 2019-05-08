# License: GNU Affero General Public License v3 or later
# A copy of GNU AGPL v3 should have been included in this software package in LICENSE.txt.

""" A collection of functions for running commands from the meme suite, FIMO and MEME.
"""

import logging

from .base import execute, get_config


def run_fimo_simple(query_motif_file: str, target_sequence: str) -> str:
    """ Runs FIMO on the provided inputs

        Arguments:
            query_motif_file: the path to the file containing query motifs
            target_sequence: the path to the file containing input sequences

        Returns:
            the output from running FIMO
    """
    command = ["fimo", "--text", "--verbosity", "1", query_motif_file, target_sequence]
    result = execute(command)
    if not result.successful():
        logging.debug('FIMO returned %d: %r while searching %r', result.return_code,
                      result.stderr, query_motif_file)
        raise RuntimeError("FIMO problem while running %s... %s" % (command, result.stderr[-100:]))
    return result.stdout


def run_fimo_version() -> str:
    """ Get the version of the fimo binary """
    fimo = get_config().executables.fimo
    command = [
        fimo,
        "-version",
    ]

    version_string = execute(command).stdout.strip()
    if not version_string:
        msg = "unexpected output from fimo: %s, check path"
        raise RuntimeError(msg % fimo)
    return version_string


def run_meme_version() -> str:
    """ Get the version of the meme binary """
    meme = get_config().executables.meme
    command = [
        meme,
        "-version",
    ]

    version_string = execute(command).stdout.strip()
    if not version_string:
        msg = "unexpected output from meme: %s, check path"
        raise RuntimeError(msg % meme)
    return version_string
