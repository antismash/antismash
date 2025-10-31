# License: GNU Affero General Public License v3 or later
# A copy of GNU AGPL v3 should have been included in this software package in LICENSE.txt.

""" Covers the alternative executable name and run-specific executable path for
    external dependencies
"""

import argparse
from collections import OrderedDict
import os
from typing import AnyStr, Dict


# since different OS package managers have different ideas on what to name things
# alternate names need to be tried
_ALTERNATE_EXECUTABLE_NAMES = {
    "diamond": [
        "diamond-aligner",  # ubuntu default
        "diamond",  # general
    ],
    "hmmpfam2": [
        "hmm2pfam",  # debian/ubuntu default
        "hmmpfam2",  # general
    ],
    "fasttree": [
        "fasttree",  # general
        "FastTree",  # as per install from source instructions
    ],
}

_NO_KNOWN_ALTS = ["hmmsearch", "hmmpress", "hmmscan", "meme", "fimo",
                  "prodigal", "blastp", "makeblastdb"]
for _binary in _NO_KNOWN_ALTS:
    assert _binary not in _ALTERNATE_EXECUTABLE_NAMES, _binary
    _ALTERNATE_EXECUTABLE_NAMES[_binary] = [_binary]


class AlternateExecutablesAction(argparse.Action):  # pylint: disable=too-few-public-methods
    """ An argparse.Action to contain multiple executables with alternate names/paths. """
    def __call__(self, parser: argparse.ArgumentParser, namespace: argparse.Namespace,  # type: ignore
                 values: AnyStr, option_string: str = None) -> None:
        if not isinstance(getattr(namespace, self.dest), argparse.Namespace):
            namespace.__dict__[self.dest] = argparse.Namespace()
        if not hasattr(namespace, self.dest):
            namespace.__dict__[self.dest] = argparse.Namespace()
        group = getattr(namespace, self.dest)
        try:
            for name, path in get_executable_paths(str(values)).items():
                setattr(group, name, path)
        except ValueError as err:
            raise argparse.ArgumentError(self, str(err))


def get_default_paths() -> Dict[str, str]:
    """ Builds a default set of paths, one for each
        executable in _ALTERNATE_EXECUTABLE_NAMES.
    """
    binaries: Dict[str, str] = OrderedDict()
    for name, alternates in _ALTERNATE_EXECUTABLE_NAMES.items():
        if name in binaries:
            continue
        path = find_executable_path(*alternates)
        if path:
            binaries[name] = path
    return binaries


def get_executable_paths(binaries_arg: str) -> Dict[str, str]:
    """ Convert a commandline option for binaries in the form:
            diamond=/full/path/to/executable,hmmpfam2=hmm2pfam
        into a Namespace, provided each binary name is in _ALTERNATE_EXECUTABLE_NAMES.
    """
    binaries: Dict[str, str] = {}
    for part in binaries_arg.split(","):
        if not part:
            continue
        subparts = part.split("=")
        if not len(subparts) == 2:
            raise ValueError(f"invalid alternate executable format: {part}")
        name, path = subparts
        if not name or not path:
            raise ValueError(f"invalid alternate executable format: {part}")
        if name not in _ALTERNATE_EXECUTABLE_NAMES:
            raise ValueError(f"unrecognised executable name: {name}")
        if os.sep in path:
            path = os.path.abspath(path)
        if os.path.isabs(path):
            if not os.path.exists(path):
                raise ValueError(f"no such file: {path!r}")
            full_path = path
        else:
            full_path = find_executable_path(path)
            if not full_path:
                raise ValueError(f"cannot find executable: {path}")
        if full_path != binaries.get(name, full_path):
            raise ValueError(f"multiple paths specified for executable: {name}")
        assert full_path
        binaries[name] = full_path
    return binaries


def find_executable_path(*names: str) -> str:
    """ Finds the full path for an executable under any of the given names
        and checks that executable permissions exist

        Arguments:
            one or more strings representing executable names/file paths to find

        Returns:
            The first of the paths that was valid or an empty string if not
    """
    if not len(names):
        raise ValueError("need to provide at least one name")
    for name in names:
        if os.path.split(name)[0] and os.path.isfile(name) and os.access(name, os.X_OK):
            return name
    for path in os.environ["PATH"].split(os.pathsep):
        for name in names:
            full_path = os.path.join(path, name)
            if os.path.isfile(full_path) and os.access(full_path, os.X_OK):
                return full_path
    return ""
