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
from .executables import get_default_paths
from .loader import load_config_from_file

_USER_FILE_NAME = os.path.expanduser('~/.antismash7.cfg')
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
            if attr in self.__dict__:
                return self.__dict__[attr]
            # for some cases of members which are their own namespace, they
            # should be created with default values if missing
            if attr == "executables":
                executables = Namespace(**get_default_paths())
                self.__dict__[attr] = executables
                return executables
            raise AttributeError(f"Config has no attribute: {attr}")

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

    def __new__(cls, namespace: Union[Namespace, Dict[str, Any]] = None) -> ConfigType:  # type: ignore
        if namespace is None:
            values: Dict[str, Any] = {}
        elif isinstance(namespace, dict):
            values = namespace
        else:
            values = namespace.__dict__
        with Config.__lock:
            if Config.__singleton is None:
                Config.__singleton = Config._Config(values)
            else:
                Config.__singleton.__dict__.update(values)
        return Config.__singleton


def update_config(values: Union[Dict[str, Any], Namespace]) -> ConfigType:
    """ Updates the Config singleton with the keys and values provided in the
        given Namespace or dict. Only intended for use in unit testing.
    """
    config = Config(values)
    assert isinstance(config, ConfigType)
    return config


def get_config(no_defaults: bool = False) -> ConfigType:
    """ Returns the current config. If the current config is empty, a default set
        will be created.

        Arguments:
            no_defaults: if True, no defaults will be created in the case of not
                         yet being set

        Returns:
            the current config
    """
    config: Union[Config, ConfigType] = Config()  # mypy doesn't like the rest of this otherwise
    assert isinstance(config, ConfigType)
    # if it's completely, build a default for easier use as a library, especially for executables
    if len(config) < 1 and not no_defaults:
        config = build_config([])
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
            with_files.append(f"@{filename}")
    with_files.extend(args)
    result = parser.parse_args(with_files)

    # if isolated, keep databases value
    if isolated:
        databases = result.database_dir
        result = parser.parse_args(args)
        result.database_dir = databases

    # set a base value for the record count limit
    default.__dict__.update({"triggered_limit": False})

    # then update with all the values from config files
    default.__dict__.update(result.__dict__)
    config = Config(default)
    assert isinstance(config, ConfigType)
    return config
