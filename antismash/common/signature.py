# License: GNU Affero General Public License v3 or later
# A copy of GNU AGPL v3 should have been included in this software package in LICENSE.txt.

""" Signature related helpers """

from typing import List

from .path import get_full_path


class Signature:
    """Secondary metabolite signature"""
    def __init__(self, name: str, _type: str, description: str, cutoff: int, path: str) -> None:
        self.name = name
        self.type = _type
        self.description = description
        self.cutoff = cutoff
        self.path = path


class HmmSignature(Signature):
    """ A container holding information on a HMM signature """
    def __init__(self, name: str, description: str, cutoff: int, hmm_path: str) -> None:
        self.hmm_file = hmm_path
        if cutoff < 0:
            raise ValueError("Signature cutoffs cannot be negative: %s" % cutoff)
        super().__init__(name, 'model', description, cutoff, self.hmm_file)


def get_signature_profiles(detail_file: str) -> List[HmmSignature]:
    """ Generates HMM signature profiles from a file.

        Paths in the file are assumed to be relative to the file itself

        Arguments:
            filename: a tab separated table, each row being a single HMM reference
                        with columns: label, description, minimum score cutoff, hmm path

        Returns:
            a list of HMMSignatures
    """
    bad_lines = []
    profiles = []
    with open(detail_file, "r") as data:
        for line in data.read().split("\n"):
            if line.startswith("#") or not line.strip():
                continue
            try:
                name, desc, cutoff, filename = line.split("\t")
            except ValueError:
                bad_lines.append(line)
                continue
            profiles.append(HmmSignature(name, desc, int(cutoff),
                                         get_full_path(detail_file, filename)))

    if bad_lines:
        raise ValueError("Invalid lines in HMM detail file (first 10):\n%s" % "\n".join(bad_lines[:10]))

    return profiles
