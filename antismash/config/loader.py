# License: GNU Affero General Public License v3 or later
# A copy of GNU AGPL v3 should have been included in this software package in LICENSE.txt.

"""Configuration handling for antiSMASH

"""
import configparser
import logging
import sys
from argparse import Namespace
from os import path

from antismash.config import update_config

_DEFAULT_NAME = 'default.cfg'
_SYS_NAME = sys.platform + '.cfg'
_USER_FILE_NAME = path.expanduser('~/.antismash.cfg')
_INSTANCE_FILE_NAME = 'instance.cfg'
_BASEDIR = path.dirname(path.abspath(__file__))


def update_config_from_file(namespace=None):
    """Load config from a default and system-specific config file and
    add it to a namespace object, but don't overwrite existing settings
    """
    if namespace is None:
        namespace = Namespace()
    default_file = path.join(_BASEDIR, _DEFAULT_NAME)
    sys_file = path.join(_BASEDIR, _SYS_NAME)
    instance_file = path.join(_BASEDIR, _INSTANCE_FILE_NAME)

    # load generic configuration settins
    config = configparser.ConfigParser()
    with open(default_file, 'r') as handle:
        config.read_file(handle)

    # load system-specific config file if available
    # also load .antismash.cfg from the user's home dir
    # and last, overriding all the other settings, instance.cfg
    config.read([sys_file, _USER_FILE_NAME, instance_file])

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

    # settings from the [DEFAULT] section go to the global namespace
    for key, value in config.items('DEFAULT'):
        key = key.replace('-', '_')
        if key not in namespace:
            namespace.__dict__[key] = value
    # store it in the singleton
    update_config(namespace)
