# License: GNU Affero General Public License v3 or later
# A copy of GNU AGPL v3 should have been included in this software package in LICENSE.txt.

"""
    Common functions for path manipulation and use
"""

import contextlib
import logging
import os
from typing import Generator, List, Optional, Union


def get_full_path(current_file: str, *args: Union[str, List[str]]) -> str:
    """ Generate the absolute path of the directory containing a file
        e.g. __file__ == os.path.join(get_full_path(__file__), filename)

        Can be given extra arguments which will be added to the result
        to generate the absolute path of a file without using both
        get_full_path and os.path.join.
        e.g. get_full_path(__file__, "data", "pfam.hmm") -> "$DIR/data/pfam.hmm"

        Arguments:
            current_file: The file from which to take the base directory from
            *args: Strings to be added to the resulting directory

        Returns:
            A string containing the fully generated path
    """
    base = os.path.dirname(os.path.abspath(current_file))
    if not args:
        return base
    extra = os.path.join(*args)  # type: ignore  # mypy can't manage the splat
    assert isinstance(extra, str)
    return os.path.join(base, extra)


def locate_file(path: str, silent: bool = False) -> Optional[str]:
    """ Checks that a given file path is valid and that read permissions exist
        for the file

        Arguments:
            path: the file path to check
            silent: if True, debug messages won't be logged

        Returns:
            The path if it was valid or None if not
    """
    if os.path.split(path)[0]:
        if os.path.isfile(path) and os.access(path, os.R_OK):
            if not silent:
                logging.debug("Found file %r", path)
            return path
    return None


@contextlib.contextmanager
def changed_directory(path: str) -> Generator:
    """ Changes working directory for the duration of the context
        e.g.:

        # in /foo/bar
        with changed_directory("/foo/baz"):
            # do some work in /foo/baz
        # automatically back to /foo/bar

        Arguments:
            path: the path of the directory to change into

        Returns:
            None
    """
    current_path = os.getcwd()
    try:
        os.chdir(path)
        yield
    finally:
        os.chdir(current_path)


def is_outdated(built_files: Union[str, List[str]], source_files: Union[str, List[str]]) -> bool:
    """ Returns True if the oldest of built_files is older than the newest of source_files

        Arguments:
            built_files: a filename, or list of filenames, of built targets
            source_files: a filename, or list of filenames, of source files used to build targets

        Returns:
            whether the oldest of built_files is older than the newest of source_files
    """
    if not built_files or not source_files:
        raise ValueError("at least one built file and one source file must be provided")

    if isinstance(built_files, str):
        built_files = [built_files]
    if isinstance(source_files, str):
        source_files = [source_files]

    for filename in built_files:
        if not locate_file(filename):
            return True

    built_time = min(os.path.getmtime(filename) for filename in built_files)
    source_time = max(os.path.getmtime(filename) for filename in source_files)
    return built_time < source_time
