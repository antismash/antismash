# License: GNU Affero General Public License v3 or later
# A copy of GNU AGPL v3 should have been included in this software package in LICENSE.txt.

"""Configuration handling for antiSMASH

"""
import configparser
import logging
import sys
from argparse import Namespace
from os import path

from antismash.config.args import Config

_config = None
_basedir = path.dirname(path.abspath(__file__))
_default_name = 'default.cfg'
_sys_name = sys.platform + '.cfg'
_user_file_name = path.expanduser('~/.antismash.cfg')
_instance_file_name = 'instance.cfg'

def update_config_from_file(namespace=None):
    """Load config from a default and system-specific config file and
    add it to a namespace object, but don't overwrite existing settings
    """
    if namespace is None:
        namespace = Namespace()
    default_file = path.join(_basedir, _default_name)
    sys_file = path.join(_basedir, _sys_name)
    instance_file = path.join(_basedir, _instance_file_name)

    # load generic configuration settins
    config = configparser.ConfigParser()
    with open(default_file, 'r') as fp:
        config.read_file(fp)

    # load system-specific config file if available
    # also load .antismash.cfg from the user's home dir
    # and last, overriding all the other settings, instance.cfg
    config.read([sys_file, _user_file_name, instance_file])

    for s in config.sections():
        if s not in namespace:
            namespace.__dict__[s] = Namespace()
        for key, value in config.items(s):
            key = key.replace('-', '_')
            if key not in namespace.__dict__[s]:
                try:
                    namespace.__dict__[s].__dict__[key] = config.getboolean(s, key)
                    continue
                except ValueError:
                    pass
                namespace.__dict__[s].__dict__[key] = value

    # settings from the [DEFAULT] section go to the global namespace
    for key, value in config.items('DEFAULT'):
        key = key.replace('-', '_')
        if key not in namespace:
            namespace.__dict__[key] = value
    # store it in the singleton
    Config(namespace)

def set_config(namespace):
    """Set a namespace object to be the global configuration"""
    logging.critical("old style set_config() used")
    Config(namespace)

def get_config():
    """Get the global configuration"""
    logging.critical("old style get_config() used")
    return Config()
