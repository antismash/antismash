# License: GNU Affero General Public License v3 or later
# A copy of GNU AGPL v3 should have been included in this software package in LICENSE.txt.

""" Helpers for type hints.
"""


class ConfigType:  # pylint: disable=too-few-public-methods
    """ Exists only to allow mypy to use the stub by this name in custom_typing.pyi """
    def __init__(self) -> None:
        raise NotImplementedError("ConfigType is a stub for typing purposes only")


class AntismashModule:  # pylint: disable=too-few-public-methods
    """ Intended only for type hints of antiSMASH modules, see custom_typing.pyi """
    def __init__(self) -> None:
        raise NotImplementedError("AntismashModule is a stub for typing purposes only")
