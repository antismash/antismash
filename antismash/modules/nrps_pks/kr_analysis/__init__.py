# License: GNU Affero General Public License v3 or later
# A copy of GNU AGPL v3 should have been included in this software package in LICENSE.txt.

""" Starcevic et al. methods for analysis of KR domains """

from antismash.common import path
from antismash.common.brawn import CacheError, ensure_alignment_cached


from .kr_analysis import KR_DOMAINS_PATH


def prepare_data(logging_only: bool = False) -> list[str]:
    """ Ensures packaged data is fully prepared

        Arguments:
            logging_only: whether to return error messages instead of raising exceptions

        Returns:
            a list of error messages (only if logging_only is True)
    """
    failures: list[str] = []
    # pre-cache alignments need to be cached and up to date
    expected_alignments = [
        KR_DOMAINS_PATH,
    ]
    cache_dir = path.get_full_path(__file__, "data")
    for fasta in expected_alignments:
        try:
            ensure_alignment_cached(fasta, cache_dir)
        except CacheError as err:
            if not logging_only:
                raise
            failures.append(str(err))
    return failures
