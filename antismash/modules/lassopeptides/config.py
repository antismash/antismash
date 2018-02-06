# License: GNU Affero General Public License v3 or later
# A copy of GNU AGPL v3 should have been included in this software package in LICENSE.txt.

"""Lassopeptide-specific configuration"""

_config = None


class LassoConfig(object):
    """Lassopeptide-specific configuration"""
    __slots__ = ('fimo_present', )

    def __init__(self):
        self.fimo_present = False


def get_config():
    """Access the lassopeptide-specific config"""
    global _config
    if _config is None:
        _config = LassoConfig()
    return _config
