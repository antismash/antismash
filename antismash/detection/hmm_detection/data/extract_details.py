#!/usr/bin/env python3
# License: GNU Affero General Public License v3 or later
# A copy of GNU AGPL v3 should have been included in this software package in LICENSE.txt.

"""
Helper script to generate a hmmdetails.txt line for the given profile
"""

from argparse import ArgumentParser, FileType
import math
import os
from typing import IO


def _main() -> None:
    parser = ArgumentParser()
    parser.add_argument("profile", type=FileType('r', encoding="utf-8"),
                        help="HMM file to extract info from")
    parser.add_argument("-n", "--name", type=str, default="",
                        help="Custom name to override the profile's name with")
    parser.add_argument("-D", "--description", type=str, default="",
                        help="Custom description to override the profile")
    parser.add_argument("-c", "--cutoff", type=int, default=-1,
                        help="Custom cutoff value, to override the profile or cover a missing TC line")
    args = parser.parse_args()

    print(run(args.profile, args.name, args.description, args.cutoff))


def run(profile: IO, name: str = "", description: str = "", cutoff: int = -1) -> str:
    """ Extracts name, description, and cutoff from an HMM profile and returns
        a string with the information that is compatible with the format used in
        'hmmdetails.txt'.

        Arguments:
            profile: an open handle to the profile content
            name: a name to use instead of any found in the profile
            description: a description to use instead of any found in the profile
            cutoff: a cutoff to use instead of any used in the profile

        Returns:
            a 'hmmdetails.txt'-compatible line of text
    """
    filename = os.path.basename(profile.name)

    line = profile.readline().strip()
    while line and line != '//' and not (name and description and (cutoff > -1)):
        line = profile.readline().strip()
        if not name and line.startswith("NAME"):
            name = line.split(" ", 1)[-1].strip()
            continue
        if not description and line.startswith("DESC"):
            description = line.split(" ", 1)[-1].strip()
            if name.startswith("TIGR") and ":" in description:
                description = description.split(":")[-1].strip()
            continue
        if cutoff < 0 and line.startswith("TC"):
            cutoff = math.floor(float(line.split(" ")[-2]))
            continue

    if cutoff < 0:
        raise RuntimeError("No TC line found and no cutoff specified")

    return "\t".join((name, description, str(cutoff), filename))


if __name__ == "__main__":
    _main()
