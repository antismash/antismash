# License: GNU Affero General Public License v3 or later
# A copy of GNU AGPL v3 should have been included in this software package in LICENSE.txt.

""" Functions and classes for managing the signatures from HMM profiles """

from typing import List

from antismash.common import path
from antismash.common import signature

_SIGNATURE_CACHE: List[signature.HmmSignature] = []


def get_signature_profiles() -> List[signature.HmmSignature]:
    """ Generates the HMM signature profiles from hmmdetails.txt
        Only does the processing once per python invocation, future runs access
        existing profiles
    """
    # if already called once, then just reuse the cached results
    if not _SIGNATURE_CACHE:
        filename = path.get_full_path(__file__, "data", "hmmdetails.txt")
        _SIGNATURE_CACHE.extend(signature.get_signature_profiles(filename))

    return list(_SIGNATURE_CACHE)
