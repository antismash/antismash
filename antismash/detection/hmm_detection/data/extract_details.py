#!/usr/bin/env python3
# License: GNU Affero General Public License v3 or later
# A copy of GNU AGPL v3 should have been included in this software package in LICENSE.txt.

from argparse import ArgumentParser, FileType
import math
import os
from typing import IO


def main() -> None:
    parser = ArgumentParser()
    parser.add_argument("profile", type=FileType('r', encoding="utf-8"), help="HMM file to extract info from")
    parser.add_argument("-n", "--name", type=str, default="", help="Custom name to override the profile's name with")
    parser.add_argument("-D", "--description", type=str, default="", help="Custom description to override the profile")
    parser.add_argument("-c", "--cutoff", type=int, default=-1, help="Custom cutoff value, to override the profile or cover a missing TC line")
    args = parser.parse_args()

    print(run(args.profile, args.name, args.description, args.cutoff))


def run(profile: IO, name: str = "", description: str = "", cutoff: int = -1) -> str:
    filename = os.path.basename(profile.name)

    line = profile.readline().strip()
    while line and line != '//' and not (name and description and (cutoff > -1) ):
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
        raise RuntimeError(f"No TC line found and no cutoff specified")

    return "\t".join((name, description, str(cutoff), filename))


if __name__ == "__main__":
    main()
