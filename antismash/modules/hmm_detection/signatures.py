# License: GNU Affero General Public License v3 or later
# A copy of GNU AGPL v3 should have been included in this software package in LICENSE.txt.

from antismash.common import path
from antismash.common.signature import Signature


class HmmSignature(Signature):
    """HMM signature"""
    def __init__(self, name, description, cutoff, hmm_filename):
        self.hmm_file = path.get_full_path(__file__, "data", hmm_filename)
        self.name = name
        if cutoff < 0:
            raise ValueError("Signature cutoffs cannot be negative: %s" % cutoff)
        super().__init__(name, 'model', description, cutoff, self.hmm_file)

def get_signature_names():
    return set(prof.name for prof in get_signature_profiles())

def get_signature_profiles():
    """ Generates the HMM signature profiles from hmmdetails.txt
        Only does the processing once per python invocation, future runs access
        existing profiles
    """
    existing = getattr(get_signature_profiles, 'existing', None)
    if existing is not None:
        return existing

    bad_lines = []

    profiles = []
    with open(path.get_full_path(__file__, "hmmdetails.txt"), "r") as data:
        for line in data.read().split("\n"):
            if line.startswith("#") or not line.strip():
                continue
            if line.count("\t") != 3:
                bad_lines.append(line)
            try:
                name, desc, cutoff, filename = line.split("\t")
            except ValueError:
                bad_lines.append(line)
                continue
            profiles.append(HmmSignature(name, desc, int(cutoff), filename))

    if bad_lines:
        raise ValueError("Invalid lines in hmmdetails:\n%s" % "\n".join(bad_lines))

    get_signature_profiles.existing = profiles

    return profiles
