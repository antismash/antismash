# License: GNU Affero General Public License v3 or later
# A copy of GNU AGPL v3 should have been included in this software package in LICENSE.txt.

""" Contains various error classes for use throughout antiSMASH """


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
