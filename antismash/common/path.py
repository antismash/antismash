# License: GNU Affero General Public License v3 or later
# A copy of GNU AGPL v3 should have been included in this software package in LICENSE.txt.

import logging
import os

def get_full_path(current_file, *args):
    "Get the full path of file_to_add in the same directory as current_file"
    base = os.path.dirname(os.path.abspath(current_file))
    extra = os.path.join(*args)
    return os.path.join(base, extra)


def locate_executable(name):
    "Find an executable in the path and return the full path"
    if os.path.split(name)[0]:
        if os.path.isfile(name) and os.access(name, os.X_OK):
            logging.debug("Found executable %r", name)
            return name
    for p in os.environ["PATH"].split(os.pathsep):
        full_name = os.path.join(p, name)
        if os.path.isfile(full_name) and os.access(full_name, os.X_OK):
            logging.debug("Found executable %r", full_name)
            return full_name

    return None


def locate_file(name):
    "Find a file and return the full path"
    if os.path.split(name)[0]:
        if os.path.isfile(name) and os.access(name, os.R_OK):
            logging.debug("Found file %r", name)
            return name
    return None
