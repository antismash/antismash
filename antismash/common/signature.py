# License: GNU Affero General Public License v3 or later
# A copy of GNU AGPL v3 should have been included in this software package in LICENSE.txt.

""" Signature related helpers """

from typing import List

from .path import get_full_path


class Signature:
    """Secondary metabolite signature"""
    def __init__(self, name: str, _type: str, description: str, cutoff: int, path: str, seed_count: int = 0) -> None:
        if name.strip() != name:
            raise ValueError(f"Signature identifiers cannot have leading or trailing whitespace: {name!r}")
        self.name = name
        self.type = _type
        self.description = description
        self.cutoff = cutoff
        self.path = path
        self.seed_count = seed_count


class HmmSignature(Signature):
    """ A container holding information on a HMM signature """
    def __init__(self, name: str, description: str, cutoff: int, hmm_path: str, seed_count: int = 0) -> None:
        self.hmm_file = hmm_path
        if cutoff < 0:
            raise ValueError(f"Signature cutoffs cannot be negative: {cutoff}")
        super().__init__(name, 'model', description, cutoff, self.hmm_file, seed_count)


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
    with open(detail_file, "r", encoding="utf-8") as data:
        for line in data.read().split("\n"):
            if line.startswith("#") or not line.strip():
                continue
            try:
                name, desc, cutoff, filename = line.split("\t")
            except ValueError:
                bad_lines.append(line)
                continue
            num_seeds = -1
            # then find and read the seeds
            with open(get_full_path(detail_file, filename), "r", encoding="utf-8") as handle:
                lines = handle.readlines()
            for line in lines:
                if line.startswith('NSEQ '):
                    num_seeds = int(line[6:].strip())
                    break
            if num_seeds == -1:
                raise ValueError(f"Unknown number of seeds for hmm file: {filename}")
            profiles.append(HmmSignature(name, desc, int(cutoff),
                                         get_full_path(detail_file, filename), num_seeds))

    if bad_lines:
        raise ValueError("Invalid lines in HMM detail file (first 10):\n%s" % "\n".join(bad_lines[:10]))

    return profiles
