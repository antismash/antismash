# License: GNU Affero General Public License v3 or later
# A copy of GNU AGPL v3 should have been included in this software package in LICENSE.txt.

""" Provides functions for interacting with PFAM databases,
    focusing on mapping of protein name to PFAM id """

import glob
import logging
import os
from typing import Dict, List

from antismash.common import path

KNOWN_MAPPINGS = {}  # type: Dict[str, Dict[str, str]] # tracks name mappings per database


def find_latest_database_version(database_dir: str) -> str:
    """ Finds the most up-to-date PFAM database version in the given directory.
        Versions are expected to be in a XY.Z format, e.g. 27.0

        Arguments:
            database_dir: a path to the antiSMASH database directory,
                    must contain a 'pfam' subdirectory with different versions

        Returns:
            the latest version number as a string, e.g. "27.0"
    """
    contents = glob.glob(os.path.join(database_dir, 'pfam', "*"))
    potentials = []
    for name in contents:
        # only names in the form 27.0, 31.0, etc are valid
        try:
            potentials.append(float(os.path.basename(name)))
        except ValueError:
            continue
    if not potentials:
        raise ValueError("No matching PFAM database in location %s" % database_dir)
    latest = sorted(potentials)[-1]
    return "%.1f" % latest


def check_db(db_path: str) -> List[str]:
    "Check that all required files exist for a database"
    failure_messages = []
    for file_name in ['Pfam-A.hmm', 'Pfam-A.hmm.h3f', 'Pfam-A.hmm.h3i',
                      'Pfam-A.hmm.h3m', 'Pfam-A.hmm.h3p']:
        if not path.locate_file(os.path.join(db_path, file_name)):
            failure_messages.append("Failed to locate file: %r in %s" % (file_name, db_path))

    return failure_messages


def get_latest_db_path(data_dir: str) -> str:
    """ Given the path to the root antiSMASH data directory, return the path
        to the latest PFAM database
    """
    return get_db_path_from_version(find_latest_database_version(data_dir), data_dir)


def get_db_version_from_path(db_path: str) -> str:
    """ Given a path to a PFAM database, pulls the version number. Paths are
        expected to match the typical antiSMASH layout of .../pfam/'VERSION'/...

        Arguments:
            db_path: the path to find

        Returns:
            the version number as a string, e.g. "31.0"
    """
    parts = db_path.split(os.sep)
    if 'pfam' not in parts:
        raise ValueError("Database path does not contain a 'pfam' directory")
    try:
        index = len(parts) - list(reversed(parts)).index("pfam") + 1
        # check the format is like a float
        _ = float(parts[index - 1])
        version = parts[index - 1]
    except ValueError:
        raise ValueError("No valid database version found in PFAM database path")

    return version


def get_db_path_from_version(version: str, data_dir: str) -> str:
    """ Given a database version and the root antiSMASH data directory, return
        the path to the specified version's main file.

        If version is "latest", the latest version will be used.
    """
    if version == "latest":
        version = find_latest_database_version(data_dir)
    return os.path.join(data_dir, "pfam", version, "Pfam-A.hmm")


def _build_mapping(database: str) -> Dict[str, str]:
    """ Build a mapping of PFAM NAME field to ACC field for a database

        Arguments:
            database: a path to the database to build a mapping of

        Returns:
            a dictionary mapping NAME field to ACC field
    """
    logging.debug("Building mapping for PFAM database version %s", get_db_version_from_path(database))
    with open(database) as handle:
        entries = handle.read().split("\n//\n")
    mapping = {}
    for entry in entries:
        if not entry:
            continue
        name, acc = [s.split()[1] for s in entry.splitlines()[1:3]]
        mapping[name] = acc
    return mapping


def get_pfam_id_from_name(name: str, database: str) -> str:
    """ Fetches a PFAM id from a protein name found in a specific database, if it exists.

        Caches mappings per database, so only the first call for a database will
        be expensive (for a PFAM database, this is approximately 3 seconds).

        Arguments:
            name: the name to find a PFAM id for
            database: the path to the database the name was found with

        Returns:
            a PFAM id as a string
    """
    if database not in KNOWN_MAPPINGS:
        KNOWN_MAPPINGS[database] = _build_mapping(database)

    return KNOWN_MAPPINGS[database][name]
