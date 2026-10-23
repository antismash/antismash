# License: GNU Affero General Public License v3 or later
# A copy of GNU AGPL v3 should have been included in this software package in LICENSE.txt.

""" Annotation-wide constants to use within record management """

import dataclasses


@dataclasses.dataclass(kw_only=True)
class _Values:
    structured_comment_key: str = "antiSMASH-Data"
    subtool_key: str = "aStool"
    tool_value: str = "antismash"


_VALUES = _Values()


def get_defaults() -> tuple[tuple[str, str], ...]:
    """ Returns the current annotation defaults """
    return tuple((f.name, getattr(_VALUES, f.name)) for f in dataclasses.fields(_VALUES))


def get_structured_comment_key() -> str:
    """ Returns the key to use in genbank structured comments.
    """
    return _VALUES.structured_comment_key


def set_structured_comment_key(key: str) -> None:
    """ Sets the key to use for genbank structured comments.
    """
    if not key:
        raise ValueError("structured comment key cannot be empty")
    _VALUES.structured_comment_key = key


def get_subtool_key() -> str:
    """ Gets the key name that will be sued when annotating a sub-tool in a
        feature qualifier, if/when the qualifier will be added for a specific
        feature.
    """
    return _VALUES.subtool_key


def set_subtool_key(name: str) -> None:
    """ Sets the key name to use when annotating a sub-tool in a feature
        qualifier, if/when the qualifier will be added for a specific feature.

        Arguments:
            name: the key to use for the qualifier
    """
    if not name:
        raise ValueError("subtool key cannot be empty")
    _VALUES.subtool_key = name


def get_tool_name() -> str:
    """ Returns the name to be annotated in the 'tool' feature qualifier,
        if/when the qualifier will be added for a specific feature.
    """
    return _VALUES.tool_value


def set_tool_name(name: str) -> None:
    """ Sets the name to be annotated in the 'tool' feature qualifier,
        if/when the qualifier will be added.

        Arguments:
            name: the value to use in the qualifier
    """
    if not name:
        raise ValueError("tool name cannot be empty")
    _VALUES.tool_value = name
