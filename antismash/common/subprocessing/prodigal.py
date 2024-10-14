# License: GNU Affero General Public License v3 or later
# A copy of GNU AGPL v3 should have been included in this software package in LICENSE.txt.

""" A collection of functions for running prodigal.
"""

from dataclasses import dataclass
import logging
from typing import Iterable

from .base import execute, get_config


@dataclass(frozen=True)
class ProdigalGene:
    """ Genes as detected by prodigal.
        Similar to locations in general, start coordinates are 0-indexed, end coordinates are 1-indexed.
    """
    name: str
    start: int
    end: int
    strand: int

    @classmethod
    def from_raw_line(cls, line: str) -> "ProdigalGene":
        """ Converts a raw prodigal output line to an instance of the class.
            e.g. '>5_7137_7874_+'
        """
        # strip the leading >
        # split the underscores
        # convert strand from +/- to int
        name, start, end, strand = line.lstrip(">").split("_")
        # convert 1-indexed start into 0-indexed for consistency throughout the whole codebase
        return cls(name, int(start) - 1, int(end), -1 if strand == "-" else 1)


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


def run_prodigal(sequence: str, options: list[str] = None) -> Iterable[ProdigalGene]:
    """ Run prodigal on the given sequence, with any additional options if provided.

        Arguments:
            sequence: the sequence to run prodigal over
            options: any additional options as provided to the command line

        Returns:
            a ProdigalGene instance for each gene found
    """
    config = get_config()

    args = [config.executables.prodigal, "-f", "sco"]
    if config.genefinding_tool == "prodigal-m" or len(sequence) < 20000:
        args.extend(["-p", "meta"])
    if options:
        args.extend(options)

    result = execute(args, stdin=f">input\n{sequence}")
    if not result.successful():
        logging.error("Failed to run prodigal: %s", result.stderr)
        raise RuntimeError(f"prodigal error: {result.stderr}")
    # remembering to remove all header lines
    return (ProdigalGene.from_raw_line(line) for line in result.stdout.splitlines()
            if line.startswith(">"))
