# License: GNU Affero General Public License v3 or later
# A copy of GNU AGPL v3 should have been included in this software package in LICENSE.txt.

""" A collection of functions for running diamond.
"""

import logging
import struct
from typing import List, Optional

from helperlibs.wrappers.io import TemporaryDirectory

from antismash.common import path

from .base import execute, get_config, RunResult


def run_diamond(subcommand: str,
                opts: Optional[List[str]] = None) -> RunResult:
    """ Run a diamond subcommand, possibly with further options.

        Arguments:
            subcommand: the diamond subcommand to run
            opts: a list of additional argument strings to pass to diamond

        Returns:
            RunResult of running diamond
    """
    config = get_config()
    if not config.executables.diamond:
        raise RuntimeError("no available diamond executable")
    with TemporaryDirectory() as temp_dir:
        params = [
            config.executables.diamond,
            subcommand,
            "--threads", str(config.cpus),
            "--tmpdir", temp_dir,
        ]

        if opts:
            params.extend(opts)

        result = execute(params)
        if not result.successful():
            message = f"diamond failed to run: {subcommand}"
            if result.stderr:
                message += f" -> {result.stderr.strip().splitlines()[-3:]}"
            raise RuntimeError(message)
    return result


def run_diamond_search(query_file: str, database_file: str, mode: str = "blastp",
                       opts: Optional[List[str]] = None) -> str:
    """ Runs diamond, comparing the given query to the given database

        Arguments:
            query_file: the path of query sequence file
            database_file: the path of the database to compare to
            mode: the mode to use (defaults to blastp)
            opts: any extra options to pass to diamond

        Returns:
            the output from running diamond
    """

    args = [
        "--db", database_file,
        "--query", query_file,
    ]

    if opts:
        args.extend(opts)

    return run_diamond(mode, args).stdout


def run_diamond_makedb(database_file: str, sequence_file: str) -> RunResult:
    """ Generate a new diamond database from a fasta file.Optional

        Arguments:
            database_file: the path to the database to generate
            sequence_file: the path to a protein FASTA file to generate the database from

        Returns:
            the RunResult running diamond
    """
    args = [
        "--db", database_file,
        "--in", sequence_file,
    ]

    return run_diamond("makedb", args)


def run_diamond_version() -> str:
    """ Get the version of the diamond executable

        Returns:
            The numeric part of "diamond version"
    """

    version_string = run_diamond("version").stdout
    if not version_string.startswith("diamond version "):
        msg = "unexpected output from diamond-executable: %s, check path"
        raise RuntimeError(msg % get_config().executables.diamond)
    # Get rid of the "diamond version" prefix
    return version_string[16:].strip()


def check_diamond_db_compatible(database_file: str) -> bool:
    """ Check if the given diamond database is compatible with the installed diamond version.

        Arguments:
            database_file: the path to the database file to check

        Returns:
            True if the database file is compatible, False otherwise
    """

    with TemporaryDirectory(change=True):
        dummy_fasta = "dummy.fa"
        dummy_db = "dummy.dmnd"
        with open(dummy_fasta, "w", encoding="utf-8") as handle:
            handle.write(">test\nM\n")
        run_diamond_makedb(dummy_db, dummy_fasta)
        compatible_format = _extract_db_format(dummy_db)

    try:
        db_format = _extract_db_format(database_file)
    except ValueError:
        return False

    if db_format != compatible_format:
        logging.debug(
            "Incompatible database format for %s. Expected %s but found %s.",
            database_file, compatible_format, db_format
        )
        return False
    return True


def _extract_db_format(database_file: str) -> int:
    """ Extract version from a diamond database file

        Arguments:
            database_file: the path to the database file to extract the version from

        Returns:
            The database version as an integer
    """
    with open(database_file, "rb") as handle:
        chunk = handle.read(16)

    if len(chunk) != 16:
        logging.debug("Database %s appears corrupted.", database_file)
        raise ValueError()

    # The uint32 in front of the format version is the build ID, we don't care
    # about it yet, but might need to care in future.

    _, db_format = struct.unpack_from("II", chunk, offset=8)

    return db_format


def check_diamond_files(definition_file: str, fasta_file: str, db_file: str,
                        logging_only: bool = False) -> List[str]:
    """ Check if the database files exist in the right version.

        Arguments:
            definition_file: the path to a database metadata file
            fasta_file: the path to a proteins fasta file
            db_file: the path to the diamond databse file
            logging_only: return a list of errors messages instead of raising errors

        Returns:
            a list of error strings
    """
    failure_messages: List[str] = []

    if path.locate_file(definition_file) is None:
        failure_messages.append(f"Failed to locate cluster definition file: {definition_file!r}")

    regen_message = ""

    if path.locate_file(fasta_file) is None:
        failure_messages.append(f"Failed to locate cluster proteins: {fasta_file!r}")
        if not logging_only:
            raise FileNotFoundError(failure_messages[-1])
    elif path.locate_file(db_file) is None:
        regen_message = f"could not find diamond database: {db_file}"
    elif not check_diamond_db_compatible(db_file):
        regen_message = f"incompatible diamond database version: {db_file}"
    elif path.is_outdated(db_file, fasta_file):
        regen_message = f"diamond database outdated: {db_file}"

    if regen_message:
        try:
            logging.debug("%s, regenerating", regen_message)
            run_diamond_makedb(db_file, fasta_file)
        except RuntimeError:
            if not logging_only:
                raise
            failure_messages.append(f"Failed to regenerate diamond database {db_file!r}")

    if failure_messages:
        failure_messages.append(f"with diamond executable: {get_config().executables.diamond}")

    return failure_messages
