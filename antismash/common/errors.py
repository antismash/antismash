# License: GNU Affero General Public License v3 or later
# A copy of GNU AGPL v3 should have been included in this software package in LICENSE.txt.

# pylint: disable=unnecessary-pass

""" Contains various error classes for use throughout antiSMASH """

from typing import Any


class AntismashError(Exception):
    """ A general catch-all for errors within antiSMASH """
    pass


class AntismashInputError(AntismashError):
    """ An error for when input to antiSMASH as a whole are invalid or contain
        unsupported content.

        Not intended to replace TypeError/ValueError in functions not directly
        parsing or converting input files.
    """
    pass


class ResultsUnusableError(AntismashInputError):
    """ An error for cases when existing results cannot be reused and cause
        unrelated modules to have unusable/misleading results.
    """
    pass


class RegionRegenerationError(ResultsUnusableError):
    """ An error explicitly for when regeneration of results would affect areas
        that were contained in previous results.
    """
    def __init__(self, msg: str, *args: Any, add_default: bool = True,
                 **kwargs: Any) -> None:
        if not msg:
            msg = "Results cannot be reused"
        if not add_default:
            super().__init__(msg, *args, **kwargs)
            return

        if msg.endswith("."):
            msg += " "
        else:
            msg += ", "
        msg += "which would affect existing regions"
        super().__init__(msg, *args, **kwargs)
