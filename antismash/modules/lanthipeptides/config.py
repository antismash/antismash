# License: GNU Affero General Public License v3 or later
# A copy of GNU AGPL v3 should have been included in this software package in LICENSE.txt.

"""Lanthipeptide-specific configuration"""


class LanthiConfig(object):  # pylint: disable=too-few-public-methods
    """Lanthipeptide-specific configuration"""
    __slots__ = ('fimo_present', )

    def __init__(self):
        self.fimo_present = False

_CONFIG = LanthiConfig()

def get_config():
    """Access the lanthieptide-specific config"""
    return _CONFIG
