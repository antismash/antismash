# License: GNU Affero General Public License v3 or later
# A copy of GNU AGPL v3 should have been included in this software package in LICENSE.txt.

import logging
import os
from typing import List

from antismash.common import path, subprocessing
from antismash.common.html_renderer import Markup, docs_link
from antismash.config import ConfigType, get_config

from .analysis import compare_precursor_cores, MultiDBResults
from .databases import get_databases

_CHECKED_DATABASES: dict[str, list[str]] = {}


def ensure_database_built(filepath: str, return_not_raise: bool = False) -> List[str]:
    """ Ensures that the given blast database exists and that the generated
        files aren't out of date.

        Arguments:
            filepath: the path to the blast database
            return_not_raise: whether to catch errors and return their messages as strings

        Returns:
            any encountered error messages, will never be populated without return_not_raise == True
    """
    components = [f"{filepath}.{ext}" for ext in ["pdb", "phr", "pin", "pot", "psq", "ptf", "pto"]]

    if path.is_outdated(components, filepath):
        logging.info(f"{filepath} components missing or obsolete, rebuilding database")
        result = subprocessing.run_makeblastdb(filepath)
        if not result.successful():
            msg = f"Failed to build blast database {filepath!r}: {result.stderr}"
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
        data_path = database.get_fasta_path(config)
        # since this repeats for many modules, just do the check once per process
        if data_path not in _CHECKED_DATABASES:
            _CHECKED_DATABASES[data_path] = ensure_database_built(data_path, logging_only)
            continue
        failure_messages.extend(_CHECKED_DATABASES[data_path])
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
            failure_messages.append(f"Failed to locate file: {binary_name!r}")

    failure_messages.extend(prepare_data(logging_only=True))

    return failure_messages


def get_tooltip_text() -> Markup:
    """ Returns generalised tooltip text, including link, for CompaRiPPson """
    return Markup(
        "Includes CompaRiPPson results for any available databases, "
        f"with a detailed explanation {docs_link('here', 'modules/comparippson')}."
    )
