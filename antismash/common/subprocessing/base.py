# License: GNU Affero General Public License v3 or later
# A copy of GNU AGPL v3 should have been included in this software package in LICENSE.txt.

""" A collection of functions and related classes for executing generic external
    commands from within antismash.
"""

import multiprocessing
import logging
import os
from subprocess import Popen, PIPE, TimeoutExpired
import sys
from typing import Any, Callable, Iterable, IO, List, Optional, Union
import warnings

from antismash.config import get_config

# Don't display the SearchIO experimental warning, we know this.
with warnings.catch_warnings():
    warnings.simplefilter("ignore")
    from Bio import SearchIO  # for import by others without wrapping, pylint: disable=unused-import


class RunResult:
    """ A container for simplifying the results of running a command """
    def __init__(self, command: List[str], stdout: bytes, stderr: bytes,
                 return_code: int, piped_out: bool, piped_err: bool) -> None:
        self.command = command
        self.stdout_piped = piped_out
        self.stderr_piped = piped_err
        if piped_out:
            self.stdout = stdout.decode()
        if piped_err:
            self.stderr = stderr.decode()
        self.return_code = return_code

    def __getattribute__(self, attr: str) -> Union[bool, str, int]:
        if attr == 'stdout' and not self.stdout_piped:
            raise ValueError("stdout was redirected to file, unable to access")
        if attr == 'stderr' and not self.stderr_piped:
            raise ValueError("stderr was redirected to file, unable to access")
        return super().__getattribute__(attr)

    def successful(self) -> bool:
        """ Returns True if the command exited with an exit code of zero """
        return not self.return_code

    def get_command_string(self) -> str:
        """ Returns the command that was run to obtain this result """
        return " ".join(self.command)


def execute(commands: List[str], stdin: Optional[str] = None, stdout: Union[int, IO[Any], None] = PIPE,
            stderr: Union[int, IO[Any], None] = PIPE, timeout: int = None) -> RunResult:
    """ Executes commands in a system-independent manner via a child process.

        By default, both stderr and stdout will be piped and the outputs
        accessible.

        Arguments:
            commands: a list of arguments to execute
            stdin: None or input to be piped into the child process
            stdout: if a file is provided, stdout from the child process
                    will be piped to that file instead of the parent process
            stderr: if a file is provided, stderr from the child process
                    will be piped to that file instead of the parent process
            timeout: if provided, the child process will be terminated after
                     this many seconds

        Returns:
            a RunResult object containing any piped output
    """
    options = get_config()
    # if config has been set up and if there is an override, use that instead
    if commands[0] in options.get("executables", {}):
        commands[0] = getattr(options.executables, commands[0])

    if stdin is not None:
        stdin_redir = PIPE  # type: Optional[int]
        input_bytes = stdin.encode("utf-8")  # type: Optional[bytes]
    else:
        stdin_redir = None
        input_bytes = None

    try:
        proc = Popen(commands, stdin=stdin_redir, stdout=stdout, stderr=stderr)
        out, err = proc.communicate(input=input_bytes, timeout=timeout)
    except TimeoutExpired:
        proc.kill()
        assert isinstance(timeout, int)
        raise RuntimeError("Child process '%s' timed out after %d seconds" % (
                commands, timeout))

    return RunResult(commands, out, err, proc.returncode, stdout == PIPE,
                     stderr == PIPE)


def parallel_function(function: Callable, args: Iterable[List[Any]],
                      cpus: Optional[int] = None, timeout: int = None) -> list:
    """ Runs the given function in parallel on `cpus` cores at a time.
        Uses separate processes so all args are effectively immutable.

        Both function and args must be picklable (i.e. the function can't be
        declared in a local scope)

        e.g. parallel_function(len, [["a"], ["aa"], ["aaa"]]) -> [1, 2, 3]

        Arguments:
            function: the function to run, cannot be a lambda
            args: a list of lists, containing the arguments for each function call
            cpus: the number of processes to start (defaults to Config.cpus)
            timeout: the maximum time to wait in seconds

        Returns:
            A list of return values of the target function.
    """

    if not cpus:
        cpus = get_config().cpus

    # if only 1 core is to be used, don't fork to run it... this ignores timeout
    if cpus == 1:
        # list() to handle generators, * to expand the list of args
        return [function(*argset) for argset in args]

    pool = multiprocessing.Pool(cpus)
    jobs = pool.starmap_async(function, args)

    timeouts = False

    try:
        results = jobs.get(timeout=timeout)
    except multiprocessing.TimeoutError:
        timeouts = True
    finally:
        pool.terminate()
        pool.join()
    if timeouts:
        raise RuntimeError("Timeout in parallel function:", function)
    return results


def child_process(command: List[str]) -> int:
    """ Called by multiprocessing's map or map_async method, cannot be locally
        defined """
    try:
        result = execute(command)
        if result.stderr:
            print(result.stderr, file=sys.stderr)
        return result.return_code
    except KeyboardInterrupt:
        #  Need to raise some runtime error that is not KeyboardInterrupt, because Python has a bug there
        raise RuntimeError("Killed by keyboard interrupt")
    return -1


def verbose_child_process(command: List[str]) -> int:
    """ A wrapper of child_process() that logs the command being run """
    logging.debug("Calling %s", " ".join(command))
    return child_process(command)


def parallel_execute(commands: List[List[str]], cpus: Optional[int] = None,
                     timeout: Optional[int] = None, verbose: bool = True) -> List[int]:
    """ Limited return vals, only returns return codes
    """
    if verbose:
        runner = verbose_child_process
    else:
        runner = child_process
    os.setpgid(0, 0)
    if not cpus:
        cpus = get_config().cpus
    assert isinstance(cpus, int)
    pool = multiprocessing.Pool(cpus)
    jobs = pool.map_async(runner, commands)

    try:
        errors = jobs.get(timeout=timeout)
    except multiprocessing.TimeoutError:
        pool.terminate()
        assert isinstance(timeout, int)
        raise RuntimeError("One of %d child processes timed out after %d seconds" % (
                cpus, timeout))

    except KeyboardInterrupt:
        logging.error("Interrupted by user")
        pool.terminate()
        raise

    pool.close()

    return errors
