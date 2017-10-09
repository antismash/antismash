# License: GNU Affero General Public License v3 or later
# A copy of GNU AGPL v3 should have been included in this software package in LICENSE.txt.

"""
    Common functions for path manipulation and use
"""

import logging
from typing import Optional
import os


def get_full_path(current_file, *args) -> str:
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
    extra = os.path.join(*args)
    return os.path.join(base, extra)


def locate_executable(name) -> Optional[str]:
    """ Find an executable in the path and return the full path

        Arguments:
            name: the name of the executable to find

        Returns:
            The absolute path to the executable as a string or None if not found
    """
    if os.path.split(name)[0]:
        if os.path.isfile(name) and os.access(name, os.X_OK):
            logging.debug("Found executable %r", name)
            return name
    for path in os.environ["PATH"].split(os.pathsep):
        full_name = os.path.join(path, name)
        if os.path.isfile(full_name) and os.access(full_name, os.X_OK):
            logging.debug("Found executable %r", full_name)
            return full_name

    return None


def locate_file(path) -> Optional[str]:
    """ Checks that a given file path is valid and that read permissions exist
        for the file

        Arguments:
            path: the file path to check

        Returns:
            The path if it was valid or None if not
    """
    if os.path.split(path)[0]:
        if os.path.isfile(path) and os.access(path, os.R_OK):
            logging.debug("Found file %r", path)
            return path
    return None
