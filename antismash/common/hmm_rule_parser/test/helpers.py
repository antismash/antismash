# License: GNU Affero General Public License v3 or later
# A copy of GNU AGPL v3 should have been included in this software package in LICENSE.txt.

# for test files, silence irrelevant and noisy pylint warnings
# pylint: disable=no-self-use,protected-access,missing-docstring,too-many-public-methods

""" helper functions for testing modules using hmm_rule_parser """

import os

from antismash.common.signature import get_signature_profiles


def check_hmm_signatures(signature_file, hmm_dir):
    sigs = get_signature_profiles(signature_file)
    for sig in sigs:
        hmm_file = os.path.abspath(os.path.join(hmm_dir, sig.path))
        assert os.path.exists(hmm_file)
        name = None
        for line in open(hmm_file):
            if line.startswith("NAME"):
                name = line.split()[-1]
        assert name
        assert name == sig.name, "%s != %s" % (name, sig.name)
