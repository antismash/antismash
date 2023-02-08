# License: GNU Affero General Public License v3 or later
# A copy of GNU AGPL v3 should have been included in this software package in LICENSE.txt.

""" Shared functions for handling brawn's cache files """

import logging
import os

from brawn import Alignment
# for simplicity of imports, also gather up a few extras
from brawn import get_aligned_pair  # pylint: disable=unused-import
from brawn.alignment import InvalidCacheFormatError, MismatchedCacheVersionError

_ALIGNMENTS = {}


class CacheError(RuntimeError):
    """ A catch-all error for reading/rebuilding brawn cache files """


def get_cached_alignment(fasta_path: str, cache_dir: str, keep_in_memory: bool = True) -> Alignment:
    """ Returns a cached brawn alignment of the given fasta

        Arguments:
            fasta_path: the path to the raw sequence alignment in FASTA format
            cache_dir: the directory containing the brawn cache file as provided
                       when the cache was built
            keep_in_memory: whether to keep the alignment in memory for future use

        Returns:
            the matching Alignment instance
    """
    real_path = os.path.join(cache_dir, f"{os.path.basename(fasta_path)}.brawn_cache")
    if real_path not in _ALIGNMENTS:
        try:
            with open(real_path, encoding="utf-8") as handle:
                alignment = Alignment.from_cache_file(handle)
        except FileNotFoundError:
            with open(fasta_path, encoding="utf-8") as handle:
                alignment = Alignment.from_file(handle)
            with open(real_path, "w", encoding="utf-8") as handle:
                alignment.to_cache_file(handle)
        if keep_in_memory:
            _ALIGNMENTS[real_path] = alignment
    else:
        alignment = _ALIGNMENTS[real_path]
        #  if it now isn't wanted, but was previously wanted, remove it
        if not keep_in_memory:
            _ALIGNMENTS.pop(real_path, None)

    return alignment


def _wrapped_unlink(path: str) -> None:
    try:
        os.unlink(path)
    except OSError as err:
        raise CacheError(f"Cannot remove stale cache file: {path}, {err}")


def ensure_alignment_cached(fasta_path: str, cache_dir: str) -> Alignment:
    """ Checks the given alignment has a valid cache file in the given directory.
        Cache files will be named the same as the original file, but with ".brawn_cache"
        appended.

        Arguments:
            fasta_path: the source alignment
            cache_dir: the directory to store the cache file

        Returns:
            the revelant brawn alignment instance

    """
    if not os.path.exists(fasta_path):
        raise CacheError(f"missing FASTA: {fasta_path}")
    cache_path = os.path.join(cache_dir, f"{os.path.basename(fasta_path)}.brawn_cache")
    # does it exist and is it valid/up to date
    if os.path.exists(cache_path):
        try:
            with open(cache_path, "r", encoding="utf-8") as handle:
                alignment = Alignment.from_cache_file(handle)
            logging.debug("Found existing cache for %s", fasta_path)
        except MismatchedCacheVersionError:
            logging.debug("Removing stale cache for %s", fasta_path)
            _wrapped_unlink(cache_path)
        except InvalidCacheFormatError:
            logging.debug("Removing invalid cache file %s", fasta_path)
            _wrapped_unlink(cache_path)

    # if it didn't exist or was removed above, rebuild it
    if not os.path.exists(cache_path):
        logging.debug("Rebuilding cache file for %s", fasta_path)
        with open(fasta_path, encoding="utf-8") as fasta:
            alignment = Alignment.from_file(fasta)
        with open(cache_path, "w", encoding="utf-8") as handle:
            alignment.to_cache_file(handle)
    return alignment
