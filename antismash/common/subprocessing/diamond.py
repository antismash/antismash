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
                opts: Optional[List[str]] = None,
                use_default_opts: bool = True) -> RunResult:
    """ Run a diamond subcommand, possibly with further options.

        Arguments:
            subcommand: the diamond subcommand to run
            opts: a list of additional argument strings to pass to diamond
            use_default_opts: use default options for diamond run (e.g. threads, tmpdir)

        Returns:
            RunResult of running diamond
    """
    config = get_config()
    if not config.executables.diamond:
        raise RuntimeError("no available diamond executable")

    def run(args: list[str]) -> RunResult:
        result = execute(args)
        if not result.successful():
            message = f"diamond failed to run: {subcommand}"
            if result.stderr:
                message += f" -> {result.stderr.strip().splitlines()[-3:]}"
            raise RuntimeError(message)
        return result

    params = [
        config.executables.diamond,
        subcommand,
    ]
    if opts:
        params.extend(opts)

    if use_default_opts:
        with TemporaryDirectory() as temp_dir:
            params.extend([
                "--threads", str(config.cpus),
                "--tmpdir", temp_dir,
            ])
            return run(params)

    return run(params)


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
    # Normally specifying a temporary directory for makedb will work, but not for some versions.
    # The bug in diamond was introduced in 2.1.0 and fixed in 2.1.7.
    # This wouldn't be a big issue, but Debian 12/bookworm has 2.1.3
    use_defaults = True
    raw_version = run_diamond_version()
    try:
        version = tuple(int(part) for part in raw_version.split("."))
    except ValueError:
        # if the version isn't in the right format, assume that it's
        # outside the versions with the bug, since those all report consistently
        version = (0, 0, 0)
    # if it's in the range of bugged versions, run without defaults while still specifying threads
    if (2, 1, 0) <= version < (2, 1, 7):
        use_defaults = False
        args.extend(["--threads", str(get_config().cpus)])

    return run_diamond("makedb", args, use_default_opts=use_defaults)


def run_diamond_version() -> str:
    """ Get the version of the diamond executable

        Returns:
            The numeric part of "diamond version"
    """

    version_string = run_diamond("version", use_default_opts=False).stdout
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
            fasta_file: the path to a fasta file of query sequences
            db_file: the path to the diamond database file
            logging_only: return a list of errors messages instead of raising errors

        Returns:
            a list of error strings
    """
    failure_messages: List[str] = []

    if path.locate_file(definition_file) is None:
        failure_messages.append(f"Failed to locate metadata file: {definition_file!r}")

    regen_message = ""

    if path.locate_file(fasta_file) is None:
        failure_messages.append(f"Failed to locate sequence file: {fasta_file!r}")
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
