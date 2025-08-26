# License: GNU Affero General Public License v3 or later
# A copy of GNU AGPL v3 should have been included in this software package in LICENSE.txt.

"""Configuration handling for antiSMASH

"""
import configparser
import os

from argparse import Namespace

from antismash.common import path

_DEFAULT_NAME = 'default.cfg'
_BASEDIR = path.get_full_path(__file__)


def load_config_from_file(default_file: str = "") -> Namespace:
    """ Load config from default config.

        Arguments:
            default_file: the path to the default config file, if not provided
                          an embedded version will be used

        Returns:
            a Namespace mapping option name to option value
    """
    namespace = Namespace()
    default_file = default_file or os.path.join(_BASEDIR, _DEFAULT_NAME)
    # load generic configuration settins
    config = configparser.ConfigParser()
    with open(default_file, "r", encoding="utf-8") as handle:
        config.read_file(handle)

    for section in config.sections():
        if section not in namespace:
            namespace.__dict__[section] = Namespace()
        for key, value in config.items(section):
            key = key.replace('-', '_')
            if key not in namespace.__dict__[section]:
                try:
                    namespace.__dict__[section].__dict__[key] = config.getboolean(section, key)
                    continue
                except ValueError:
                    pass
                namespace.__dict__[section].__dict__[key] = value
    if "top" in namespace:
        top_level = namespace.__dict__.pop("top")
        for key, value in top_level.__dict__.items():
            namespace.__dict__[key] = value
    else:
        namespace.__dict__["branding"] = "antiSMASH"

    # settings from the [DEFAULT] section go to the global namespace
    for key, value in config.items('DEFAULT'):
        key = key.replace('-', '_')
        if key not in namespace:
            namespace.__dict__[key] = value

    return namespace
