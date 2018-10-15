# License: GNU Affero General Public License v3 or later
# A copy of GNU AGPL v3 should have been included in this software package in LICENSE.txt.

""" Logging related code too involved to inline into other parts of the code. """

import contextlib
import logging
import os
import sys
from typing import Any, Generator
from typing import Dict  # comment hints, pylint: disable=unused-import


@contextlib.contextmanager
def changed_logging(logfile: str = None, verbose: bool = False, debug: bool = False) -> Generator:
    """ Changes logging setup for the duration of the context
        e.g.:

        with changed_logging(logfile="test.log", debug=True):
            logging.warning("warning")  # will appear on console and in test.log
            logging.debug("debug")  # will appear on console and in test.log
        logging.warning("warning")  # will appear on console only
        logging.debug("debug")  # will not appear in either

        Arguments:
            logfile: None or the path to a file to write logging messages to
            verbose: whether to show INFO level messages and above
            debug: whether to show DEBUG level messages and above

        Returns:
            None
    """
    try:
        def new_critical(*args: Any) -> None:
            """ make critical messages yellow and without the normal timestamp """
            msg = "\033[1;33m{}\033[0m".format(args[0])
            print(msg % args[1:], file=sys.stderr)
        logging.critical = new_critical  # type: ignore

        log_level = logging.WARNING
        if debug:
            log_level = logging.DEBUG
        elif verbose:
            log_level = logging.INFO

        log_format = '%(levelname)-8s %(asctime)s   %(message)s'
        date_format = "%d/%m %H:%M:%S"
        logging.basicConfig(format=log_format, datefmt=date_format)
        logger = logging.getLogger()
        original_log_level = logger.getEffectiveLevel()
        logger.setLevel(log_level)

        handler = None
        original_levels = {}  # type: Dict[logging.Handler, int]

        if logfile:
            # since INFO is always wanted in the logfile, set the global level to that
            if log_level > logging.INFO:
                logger.setLevel(logging.INFO)
                # and make sure the existing handlers are more restrictive
                for stream in logger.handlers:
                    original_levels[stream] = stream.level
                    stream.setLevel(log_level)

            # create any container directory if required
            dirname = os.path.dirname(logfile)
            if dirname and not os.path.exists(dirname):
                os.makedirs(dirname)
            handler = logging.FileHandler(logfile)
            # always treat the logfile as at least verbose
            handler.setLevel(min(log_level, logging.INFO))
            handler.setFormatter(logging.Formatter(fmt=log_format, datefmt=date_format))
            logger.addHandler(handler)
        yield
    finally:
        # put things back into original condition
        logger.setLevel(original_log_level)
        if handler:
            # remove the file handler and ensure it's closed
            logger.removeHandler(handler)
            handler.flush()
            handler.close()
            # and set all other handlers back to their original level
            for stream, original_level in original_levels.items():
                stream.setLevel(original_level)
