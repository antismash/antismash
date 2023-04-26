# License: GNU Affero General Public License v3 or later
# A copy of GNU AGPL v3 should have been included in this software package in LICENSE.txt.

""" A collection of functions for running commands from the meme suite, FIMO and MEME.
"""

from dataclasses import dataclass, fields
import logging
from typing import Union

from .base import execute, get_config


@dataclass
class FIMOMotif:
    """ A version-ambiguous container for FIMO results """
    pattern_name: str
    alternate_name: str
    sequence_name: str
    start: int
    stop: int
    strand: str
    score: float
    p_value: float
    q_value: float
    matched_sequence: str

    @classmethod
    def _from_parts(cls, parts: list[str]) -> "FIMOMotif":
        """ Converts the line chunks to the relevant type for the particular field,
            then constructs the class """
        kwargs: dict[str, Union[str, int, float]] = {}
        for part, field in zip(parts, fields(cls)):
            if field.type in (int, float):
                part = field.type(part or 0)
            kwargs[field.name] = part
        return cls(**kwargs)  # type: ignore  # this is validated above and mypy can't handle it

    @classmethod
    def from_legacy_line(cls, line: str) -> "FIMOMotif":
        """ Converts a line from legacy output (< 4.11.3) to an instance """
        parts = line.strip("\n").split("\t")
        parts.insert(1, "")
        return cls._from_parts(parts)

    @classmethod
    def from_output_line(cls, line: str) -> "FIMOMotif":
        """ Converts a line from output to an instance """
        parts = line.strip("\n").split("\t")
        assert parts[0][0] != "#" and parts[0] != "motif_id", line
        if len(parts) < 10:  # tne it's legacy format
            return cls.from_legacy_line(line)
        return cls._from_parts(parts)


def read_fimo_output(output: str) -> list[FIMOMotif]:
    """ Reads FIMO output without needing to know which version of FIMO was used
        to generate the output.

        Arguments:
            output: the FIMO output to parse

        Returns:
            a list of FIMOMotif instances
    """
    motifs = []
    for line in output.splitlines():
        if line.startswith("#") or line.startswith("motif_id") or not line:
            continue
        motifs.append(FIMOMotif.from_output_line(line))
    return motifs


def run_fimo_simple(query_motif_file: str, target_sequence: str) -> list[FIMOMotif]:
    """ Runs FIMO on the provided inputs

        Arguments:
            query_motif_file: the path to the file containing query motifs
            target_sequence: the path to the file containing input sequences

        Returns:
            a list of motifs found by FIMO
    """
    command = ["fimo", "--text", "--verbosity", "1", query_motif_file, target_sequence]
    result = execute(command)
    if not result.successful():
        logging.debug('FIMO returned %d: %r while searching %r', result.return_code,
                      result.stderr, query_motif_file)
        raise RuntimeError(f"FIMO problem while running {command}... {result.stderr[-100:]}")
    return read_fimo_output(result.stdout)


def run_fimo_version() -> str:
    """ Get the version of the fimo binary """
    fimo = get_config().executables.fimo
    command = [
        fimo,
        "-version",
    ]

    version_string = execute(command).stdout.strip()
    if not version_string:
        raise RuntimeError(f"unexpected output from fimo: {fimo}, check path")
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
