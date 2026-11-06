#!/usr/bin/env python3
# License: GNU Affero General Public License v3 or later
# A copy of GNU AGPL v3 should have been included in this software package in LICENSE.txt.

""" Helper script to generate a hmmdetails.txt line for the given subcluster profile

Output columns (tab-separated): name  description  cutoff  hmm_file  [accession]
"""

from argparse import ArgumentParser, FileType
import math
import os
from typing import IO


def _main() -> None:
    parser = ArgumentParser(description=__doc__)
    parser.add_argument("profile", type=FileType("r", encoding="utf-8"),
                        help="HMM file to extract info from")
    parser.add_argument("-n", "--name", type=str, default="",
                        help="Custom name to override the profile's NAME field")
    parser.add_argument("-D", "--description", type=str, default="",
                        help="Custom description to override the profile's DESC field")
    parser.add_argument("-c", "--cutoff", type=int, default=-1,
                        help="Custom cutoff, overrides TC line or covers a missing one")
    parser.add_argument("-a", "--accession", type=str, default="",
                        help="Pfam accession to append as 5th column (e.g. PF00155)")
    args = parser.parse_args()

    print(run(args.profile, args.name, args.description, args.cutoff, args.accession))


def run(profile: IO, name: str = "", description: str = "", cutoff: int = -1,
        accession: str = "") -> str:
    """ Extracts name, description, cutoff, and optionally Pfam accession from
        an HMM profile and returns a tab-separated line compatible with the
        subcluster hmmdetails.txt format.

        Arguments:
            profile: an open handle to the profile content
            name: a name to use instead of any found in the profile
            description: a description to use instead of any found in the profile
            cutoff: a cutoff to use instead of any found in the profile
            accession: Pfam accession override; if empty, read from ACC field

        Returns:
            a hmmdetails.txt-compatible line of text
    """
    # profiles are kept in the rule_hmms directory next to the details file
    filename = f"rule_hmms/{os.path.basename(profile.name)}"
    found_acc = ""

    line = profile.readline().strip()
    while line and line != "//":
        if not name and line.startswith("NAME"):
            name = line.split(None, 1)[-1].strip()
        elif not found_acc and line.startswith("ACC"):
            found_acc = line.split(None, 1)[-1].strip()
        elif not description and line.startswith("DESC"):
            description = line.split(None, 1)[-1].strip()
            if name.startswith("TIGR") and ":" in description:
                description = description.split(":")[-1].strip()
        elif cutoff < 0 and line.startswith("TC"):
            cutoff = math.floor(float(line.split()[-2]))
        if name and description and cutoff >= 0 and found_acc:
            break
        line = profile.readline().strip()

    if not accession:
        accession = found_acc

    if cutoff < 0:
        raise RuntimeError("No TC line found and no cutoff specified")

    parts = [name, description, str(cutoff), filename]
    if accession:
        parts.append(accession)
    return "\t".join(parts)


if __name__ == "__main__":
    _main()
