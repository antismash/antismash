# License: GNU Affero General Public License v3 or later
# A copy of GNU AGPL v3 should have been included in this software package in LICENSE.txt.

""" Functions and classes for managing the signatures from HMM profiles """

from typing import List

from antismash.common import path
from antismash.common import signature


def get_signature_profiles() -> List[signature.HmmSignature]:
    """ Generates the HMM signature profiles from hmmdetails.txt
        Only does the processing once per python invocation, future runs access
        existing profiles
    """
    # if already called once, then just reuse the cached results
    existing = getattr(get_signature_profiles, 'existing', None)
    if existing is not None:
        assert isinstance(existing, list)
        return existing

    # not cached yet, so generate
    profiles = signature.get_signature_profiles(path.get_full_path(__file__, "data", "hmmdetails.txt"))

    # cache this for future reuse, and silence mypy warnings because it can't handle it
    get_signature_profiles.existing = profiles  # type: ignore

    return profiles
