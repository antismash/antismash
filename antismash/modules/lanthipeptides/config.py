# License: GNU Affero General Public License v3 or later
# A copy of GNU AGPL v3 should have been included in this software package in LICENSE.txt.

"""Lanthipeptide-specific configuration"""

_config = None


class LanthiConfig(object):
    """Lanthipeptide-specific configuration"""
    __slots__ = ('fimo_present', )

    def __init__(self):
        self.fimo_present = False


def get_config():
    """Access the lanthieptide-specific config"""
    global _config
    if _config is None:
        _config = LanthiConfig()
    return _config
