# License: GNU Affero General Public License v3 or later
# A copy of GNU AGPL v3 should have been included in this software package in LICENSE.txt.

import threading

def update_config(values):
    return Config(values)

def get_config():
    return Config()

def destroy_config() -> None:
    Config().__dict__.clear()


class Config:  # since it's a glorified namespace, pylint: disable=too-few-public-methods
    __singleton = None
    __lock = threading.Lock()
    class _Config():
        def __init__(self, indict):
            if indict:
                self.__dict__.update(indict)

        def get(self, key, default=None):
            return self.__dict__.get(key, default)

        def __getattr__(self, attr):
            return self.__dict__[attr]

        def __setattr__(self, attr, value):
            raise RuntimeError("Config options can't be set directly")

        def __iter__(self):
            for i in self.__dict__.items():
                yield i

        def __repr__(self):
            return str(self)

        def __str__(self):
            return str(dict(self))

        def __len__(self):
            return len(self.__dict__)

    def __new__(cls, namespace=None):
        if namespace is None:
            values = {}
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
