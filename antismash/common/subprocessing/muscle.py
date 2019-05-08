# License: GNU Affero General Public License v3 or later
# A copy of GNU AGPL v3 should have been included in this software package in LICENSE.txt.

""" A collection of functions for running muscle.
"""

from tempfile import NamedTemporaryFile
from typing import Dict

from antismash.common.fasta import read_fasta, write_fasta

from .base import execute, get_config


def run_muscle_single(seq_name: str, seq: str, comparison_file: str) -> Dict[str, str]:
    """ Runs muscle over a single sequence against a comparison file in profile
        mode and returns a dictionary of the resulting alignments

        Arguments:
            seq_name: the name of the query
            seq: the sequence to align
            comparison_file: the path of the file containing comparison sequences

        Returns:
            a dictionary mapping sequence name (query or reference) to alignment
    """
    with NamedTemporaryFile(mode="w+") as temp_in:
        with NamedTemporaryFile(mode="w+") as temp_out:
            write_fasta([seq_name], [seq], temp_in.name)
            # Run muscle and collect sequence positions from file
            result = execute([get_config().executables.muscle,
                              "-profile", "-quiet",
                              "-in1", comparison_file,
                              "-in2", temp_in.name,
                              "-out", temp_out.name])
            if not result.successful():
                raise RuntimeError("muscle returned %d: %r while comparing query named %s" % (
                                   result.return_code, result.stderr.replace("\n", ""),
                                   seq_name))
            fasta = read_fasta(temp_out.name)
    return fasta


def run_muscle_version() -> str:
    """ Get the version of the muscle binary """
    muscle = get_config().executables.muscle
    command = [
        muscle,
        "-version",
    ]

    version_string = execute(command).stdout
    if not version_string.startswith("MUSCLE"):
        msg = "unexpected output from muscle: %s, check path"
        raise RuntimeError(msg % muscle)
    # get rid of the non-version stuff in the output
    return version_string.split()[1]
