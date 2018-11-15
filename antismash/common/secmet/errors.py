# License: GNU Affero General Public License v3 or later
# A copy of GNU AGPL v3 should have been included in this software package in LICENSE.txt.

""" Contains Exception subclasses for more specific communication of issues """


class SecmetInvalidInputError(ValueError):
    """ An error for when input is invalid or contains unsupported content.

        Not intended to replace TypeError/ValueError in functions not directly
        parsing or converting input files.
    """
    pass
