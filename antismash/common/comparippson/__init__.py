# License: GNU Affero General Public License v3 or later
# A copy of GNU AGPL v3 should have been included in this software package in LICENSE.txt.

import logging
import os
from typing import List

from antismash.common import path, subprocessing
from antismash.config import ConfigType, get_config

from .analysis import compare_precursor_cores, MultiDBResults
from .databases import get_databases


def ensure_database_built(filepath: str, return_not_raise: bool = False) -> List[str]:
    """ Ensures that the given blast database exists and that the generated
        files aren't out of date.

        Arguments:
            filepath: the path to the blast database
            return_not_raise: whether to catch errors and return their messages as strings

        Returns:
            any encountered error messages, will never be populated without return_not_raise == True
    """
    components = ["{}.{}".format(filepath, ext) for ext in ["pdb", "phr", "pin", "pot", "psq", "ptf", "pto"]]

    if path.is_outdated(components, filepath):
        logging.info(f"{filepath} components missing or obsolete, rebuilding database")
        result = subprocessing.run_makeblastdb(filepath)
        if not result.successful():
            msg = "Failed to build blast database {!r}: {}".format(filepath, result.stderr)
            if not return_not_raise:
                raise RuntimeError(msg)
            return [msg]
    return []


def prepare_data(logging_only: bool = False) -> List[str]:
    """ Prepare the relevant databases. Expected directory layout is:
         'comparippson/category/version/cores.fa'

        Arguments:
            logging_only: whether to log errors or raise errors

        Returns:
            a list of error messages, one for each failure (if any)
    """
    failure_messages: List[str] = []
    config = get_config()

    data_dir = os.path.join(config.database_dir, "comparippson")
    if "mounted_at_runtime" in data_dir:  # can't prepare these
        return failure_messages

    for database in get_databases(config):
        failure_messages.extend(ensure_database_built(database.get_fasta_path(config), return_not_raise=logging_only))

    return failure_messages


def check_prereqs(options: ConfigType) -> List[str]:
    """ Checks if all required binaries are available and that all databases are built

        Arguments:
            options: antiSMASH config object

        Returns:
            a list of error messages, one for each failure (if any)
    """
    required_binaries = [
        'blastp',
    ]

    failure_messages = []
    for binary_name in required_binaries:
        if binary_name not in options.executables:
            failure_messages.append("Failed to locate file: %r" % binary_name)

    failure_messages.extend(prepare_data(logging_only=True))

    return failure_messages
