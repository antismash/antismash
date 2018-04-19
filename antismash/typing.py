# License: GNU Affero General Public License v3 or later
# A copy of GNU AGPL v3 should have been included in this software package in LICENSE.txt.

""" Helpers for type hints.
"""


class AntismashModule:  # pylint: disable=too-few-public-methods
    """ Intended only for type hints of antiSMASH modules, see typing.pyi """
    def __init__(self) -> None:
        raise NotImplementedError("Antismash module is a stub for typing purposes only")
