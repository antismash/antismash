# License: GNU Affero General Public License v3 or later
# A copy of GNU AGPL v3 should have been included in this software package in LICENSE.txt.

""" Manages commandline/runtime options for antiSMASH

    Options are available at any given moment by use of get_config(), this makes
    use of a singleton object.

    For the purposes of testing only, update_config() and destroy_config() are
    provided.

"""


from argparse import Namespace
import os
import threading
from typing import Any, Dict, Iterator, List, Optional, Tuple, Union

from antismash.custom_typing import AntismashModule, ConfigType

from .args import build_parser, AntismashParser
from .loader import load_config_from_file

_USER_FILE_NAME = os.path.expanduser('~/.antismash5.cfg')
_INSTANCE_FILE_NAME = os.path.abspath(os.path.join(os.path.dirname(__file__), 'instance.cfg'))


class Config:  # since it's a glorified namespace, pylint: disable=too-few-public-methods
    """ Keeps options values for antiSMASH.

        Really just the constructor for an internal singleton keeping state
    """
    __singleton = None
    __lock = threading.Lock()

    class _Config(ConfigType):  # ConfigType used for typing stub only
        """ The real options object.
            Only one of these should exist for the lifetime of an antiSMASH run.
        """
        def __init__(self, indict: Dict[str, Any]) -> None:  # pylint: disable=super-init-not-called
            if indict:
                self.__dict__.update(indict)

        def get(self, key: str, default: Any = None) -> Any:
            """ Returns a value for a key if the key exists, otherwise returns
                the default value provided or None """
            return self.__dict__.get(key, default)

        def __getattr__(self, attr: str) -> Any:
            return self.__dict__[attr]

        def __setattr__(self, attr: str, value: Any) -> None:
            # special exceptions for those arguments we must update after
            # reading sequences and before analysis
            if attr in ["output_dir", "version", "all_enabled_modules"]:
                self.__dict__[attr] = value
                return
            raise RuntimeError("Config options can't be set directly")

        def __iter__(self) -> Iterator[Tuple[str, Any]]:
            for i in self.__dict__.items():
                yield i

        def __repr__(self) -> str:
            return str(self)

        def __str__(self) -> str:
            return str(dict(self))

        def __len__(self) -> int:
            return len(self.__dict__)

    def __new__(cls, namespace: Union[Namespace, Dict[str, Any]] = None) -> ConfigType:
        if namespace is None:
            values = {}  # type: Dict[str, Any]
        elif isinstance(namespace, dict):
            values = namespace
        else:
            values = namespace.__dict__
        Config.__lock.acquire()
        if Config.__singleton is None:
            Config.__singleton = Config._Config(values)
        else:
            Config.__singleton.__dict__.update(values)
        Config.__lock.release()
        return Config.__singleton


def update_config(values: Union[Dict[str, Any], Namespace]) -> ConfigType:
    """ Updates the Config singleton with the keys and values provided in the
        given Namespace or dict. Only intended for use in unit testing.
    """
    config = Config(values)
    assert isinstance(config, ConfigType)
    return config


def get_config() -> ConfigType:
    """ Returns the current config """
    config = Config()
    assert isinstance(config, ConfigType)
    return config


def destroy_config() -> None:
    """ Destroys all settings in the Config singleton. Only intended for use
        with unit testing and when antismash is loaded as a library for multiple
        runs with differing options.
    """
    Config().__dict__.clear()


def build_config(args: List[str], parser: Optional[AntismashParser] = None, isolated: bool = False,
                 modules: List[AntismashModule] = None) -> ConfigType:
    """ Builds up a Config. Uses, in order of lowest priority, a users config
        file (~/.antismash5.cfg), an instance config file
        (antismash/config/instance.cfg), and the provided command line options.

        If in isolated mode, only database directory
    even if isolated, the database directory is vital, so load the files
        and keep the database value, unless it will be overridden on the command
        line
    """
    # load default for static information, e.g. URLs
    default = load_config_from_file()

    if not parser:
        parser = build_parser(from_config_file=True, modules=modules)
    # check if the simple version will work
    if isolated and "--databases" in args:
        default.__dict__.update(parser.parse_args(args).__dict__)
        config = Config(default)
        assert isinstance(config, ConfigType)
        return config

    # load from file regardless because we need databases
    with_files = []
    for filename in (_USER_FILE_NAME, _INSTANCE_FILE_NAME):
        if os.path.exists(filename):
            with_files.append("@%s" % filename)
    with_files.extend(args)
    result = parser.parse_args(with_files)

    # if isolated, keep databases value
    if isolated:
        databases = result.database_dir
        result = parser.parse_args(args)
        result.database_dir = databases

    default.__dict__.update(result.__dict__)
    config = Config(default)
    assert isinstance(config, ConfigType)
    return config
