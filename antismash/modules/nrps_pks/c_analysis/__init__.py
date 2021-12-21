# License: GNU Affero General Public License v3 or later
# A copy of GNU AGPL v3 should have been included in this software package in LICENSE.txt.

""" Analysis of Condensation domains """

import os

from antismash.common import path
from antismash.common.brawn import CacheError, ensure_alignment_cached


from .c_analysis import DATA_DIR, C_DOMAINS_FILENAME


def prepare_data(logging_only: bool = False) -> list[str]:
    """ Ensures packaged data is fully prepared

        Arguments:
            logging_only: whether to return error messages instead of raising exceptions

        Returns:
            a list of error messages (only if logging_only is True)
    """
    failures: list[str] = []
    try:
        ensure_alignment_cached(os.path.join(DATA_DIR, C_DOMAINS_FILENAME), DATA_DIR)
    except CacheError as err:
        if not logging_only:
            raise
        failures.append(str(err))
    return failures
