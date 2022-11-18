# License: GNU Affero General Public License v3 or later
# A copy of GNU AGPL v3 should have been included in this software package in LICENSE.txt.

""" Provides functions for interacting with PFAM databases,
    focusing on mapping of protein name to PFAM id """

import logging
import os
from typing import Dict, List, Optional

from antismash.common import hmmer, path

KNOWN_MAPPINGS: Dict[str, Dict[str, str]] = {}  # tracks name mappings per database
KNOWN_CUTOFFS: Dict[str, Dict[str, float]] = {}  # tracks profile cutoffs per database


def find_latest_database_version(database_dir: str) -> str:
    """ Finds the most up-to-date PFAM database version in the given directory.
        Versions are expected to be in a XY.Z format, e.g. 27.0

        Arguments:
            database_dir: a path to the antiSMASH database directory,
                    must contain a 'pfam' subdirectory with different versions

        Returns:
            the latest version number as a string, e.g. "27.0"
    """
    return path.find_latest_database_version(os.path.join(database_dir, "pfam"),
                                             required_file_pattern="Pfam-A.hmm")


def check_db(db_path: str) -> List[str]:
    "Check that all required files exist for a database"
    failure_messages = []
    file_name = 'Pfam-A.hmm'
    full_path = os.path.join(db_path, file_name)
    if not path.locate_file(full_path):
        failure_messages.append("Failed to locate file: %r in %s" % (file_name, db_path))
    else:
        failure_messages.extend(hmmer.ensure_database_pressed(full_path, return_not_raise=True))
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
    logging.debug("Building name mapping for HMMer database %s", os.path.basename(database))

    with open(database) as handle:
        entries = handle.read().split("\n//\n")
    mapping = {}
    for entry in entries:
        if not entry:
            continue
        name, acc = [s.split()[1] for s in entry.splitlines()[1:3]]
        mapping[name] = acc
    return mapping


def _build_cutoff_mapping(database: str) -> Dict[str, float]:
    """ Build a mapping of Pfam ACC field to TC field for a database

        Arguments:
            database: a path to the database to build a mapping from

        Returns:
            a dictionary mapping ACC field to TC value
    """
    logging.debug("Building cutoff mapping for HMMer database %s", database)
    with open(database) as handle:
        entries = handle.read().split("\n//\n")
    mapping = {}
    for entry in entries:
        if not entry:
            continue
        lines = iter(entry.splitlines())
        acc: Optional[str] = None
        cutoff: Optional[float] = None
        try:
            while acc is None:
                line = next(lines)
                if line.startswith("ACC "):
                    try:
                        acc = line.split()[1]
                    except IndexError:
                        raise ValueError(f"profile accession line malformed: {line}")
                    break
            while cutoff is None:
                line = next(lines)
                if line.startswith("TC "):
                    try:
                        cutoff = float(line.split()[1])
                    except IndexError:
                        raise ValueError(f"profile cutoff line malformed {acc}: {line}")
                    except TypeError as err:
                        raise TypeError(f"profile cutoff invalid for {acc}: {line}") from err
                    break
        except StopIteration:
            if acc:
                raise ValueError(f"profile missing threshold cutoff: {acc}")
            raise ValueError("profile missing accession")
        assert acc and cutoff
        mapping[acc] = cutoff
    return mapping


def get_pfam_cutoffs(database: str) -> Dict[str, float]:
    """ Returns a mapping of accession to trusted cutoff for the given database.

        Caches results for performance.

        Arguments:
            database: the path to the database to gather cutoffs from

        Returns:
            a dictionary mapping ACC field to TC value
    """
    if database not in KNOWN_CUTOFFS:
        KNOWN_CUTOFFS[database] = _build_cutoff_mapping(database)

    return KNOWN_CUTOFFS[database]


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
