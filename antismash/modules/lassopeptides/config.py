# License: GNU Affero General Public License v3 or later
# A copy of GNU AGPL v3 should have been included in this software package in LICENSE.txt.

"""Lassopeptide-specific configuration"""


class LassoConfig:  # pylint: disable=too-few-public-methods
    """Lassopeptide-specific configuration"""
    __slots__ = ('fimo_present', )

    def __init__(self) -> None:
        self.fimo_present = False


_CONFIG = LassoConfig()


def get_config() -> LassoConfig:
    """Access the lassopeptide-specific config"""
    return _CONFIG
