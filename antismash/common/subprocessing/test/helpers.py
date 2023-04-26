# License: GNU Affero General Public License v3 or later
# A copy of GNU AGPL v3 should have been included in this software package in LICENSE.txt.

# for test files, silence irrelevant and noisy pylint warnings
# pylint: disable=use-implicit-booleaness-not-comparison,protected-access,missing-docstring

from antismash.common import subprocessing


class DummyConfig:
    """ A dummy object to mimc the Config object, complete with executables """
    class DummyExecutables:
        def __init__(self, paths: dict[str, str]):
            for key, value in paths.items():
                setattr(self, key, value)

    def __init__(self, cpus=1, paths=None):
        paths = paths or {}
        self.executables = DummyConfig.DummyExecutables(paths)
        self.cpus = cpus


class DummyResult(subprocessing.RunResult):
    """ A dummy object for mocking subprocessing.base.execute """
    def __init__(self, stdout: str = "", stderr: str = "", return_code: int = 0):
        super().__init__(["dummy"], stdout.encode(), stderr.encode(),
                         return_code=return_code, piped_out=True, piped_err=True)
