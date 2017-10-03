# License: GNU Affero General Public License v3 or later
# A copy of GNU AGPL v3 should have been included in this software package in LICENSE.txt.

from io import StringIO
import logging
import multiprocessing
import os
from subprocess import PIPE, Popen, TimeoutExpired

import warnings
# Don't display the SearchIO experimental warning, we know this.
with warnings.catch_warnings():
    warnings.simplefilter("ignore")
    from Bio import SearchIO

from helperlibs.wrappers.io import TemporaryDirectory

from antismash.config import get_config

class RunResult:
    def __init__(self, command, stdout, stderr, return_code, piped_out, piped_err):
        self.command = command
        self.stdout_piped = piped_out
        self.stderr_piped = piped_err
        if piped_out:
            self.stdout = stdout.decode()
        if piped_err:
            self.stderr = stderr.decode()
        self.return_code = return_code

    def __getattribute__(self, attr):
        if attr == 'stdout' and not self.stdout_piped:
            raise ValueError("stdout was redirected to file, unable to access")
        if attr == 'stderr' and not self.stderr_piped:
            raise ValueError("stderr was redirected to file, unable to access")
        return super().__getattribute__(attr)

    def successful(self) -> int:
        return not self.return_code

    def get_command_string(self) -> str:
        return " ".join(self.command)

def execute(commands, stdin=None, stdout=PIPE, stderr=PIPE, timeout=None) -> RunResult:
    "Execute commands in a system-independent manner"

    if stdin is not None:
        stdin_redir = PIPE
        input_bytes = stdin.encode("utf-8")
    else:
        stdin_redir = None
        input_bytes = None

    try:
        proc = Popen(commands, stdin=stdin_redir, stdout=stdout, stderr=stderr)
        out, err = proc.communicate(input=input_bytes, timeout=timeout)
    except TimeoutExpired:
        proc.kill()
        raise RuntimeError("Child process '%s' timed out after %d seconds" % (
                commands, timeout))

    return RunResult(commands, out, err, proc.returncode, stdout == PIPE,
                     stderr == PIPE)

def _helper(args):
    """ Needed for parallel_function to split function to call and args.

        Since functions passed to parallel_function can't be locally defined,
        this has to be declared here
    """
    return args[0](*args[1:])

def parallel_function(function, args, cpus=None, timeout=None) -> list:
    """ Runs the given function in parallel on `cpus` cores at a time.
        Uses separate processes so all args are effectively immutable.

        Both function and args must be picklable

        Arguments:
            function: the function to run, cannot be a lambda
            args: the arguments for each function call
            cpus: the number of processes to start (defaults to Config.cpus)
            timeout: the maximum time to wait in seconds

        Returns:
            A list of return values of the target function.

        e.g. pool_parallel_function(len, ["a", "aa", "aaa"]) -> [1, 2, 3]
    """

    if not cpus:
        cpus = get_config().cpus
    p = multiprocessing.Pool(cpus)
    jobs = p.map_async(_helper, ([function] + arglist for arglist in args))

    timeouts = False

    try:
        results = jobs.get(timeout=timeout)
    except multiprocessing.TimeoutError:
        timeouts = True
    finally:
        p.terminate()
        p.join()
    if timeouts:
        raise RuntimeError("Timeout in parallel function:", function)
    return results


def child_process(command):
    """ Called by multiprocessing's map or map_async method, cannot be locally
        defined """
    try:
        logging.debug("Calling %s", " ".join(command))
        return execute(command).return_code
    except KeyboardInterrupt:
        #  Need to raise some runtime error that is not KeyboardInterrupt, because Python has a bug there
        raise RuntimeError("Killed by keyboard interrupt")
    return -1

def parallel_execute(commands, cpus=None, timeout=None):
    """ Limited return vals, only returns return codes
    """
    os.setpgid(0, 0)
    if not cpus:
        cpus = get_config().cpus
    p = multiprocessing.Pool(cpus)
    jobs = p.map_async(child_process, commands)

    try:
        errors = jobs.get(timeout=timeout)
    except multiprocessing.TimeoutError:
        p.terminate()
        raise RuntimeError("One of %d child processes timed out after %d seconds" % (
                cpus, timeout))
    except KeyboardInterrupt:
        logging.error("Interrupted by user")
        p.terminate()
        raise

    p.close()

    return errors

def run_hmmsearch(query_hmmfile, target_sequence, use_tempfile=False):
    "Run hmmsearch"
    config = get_config()
    command = ["hmmsearch", "--cpu", str(config.cpus),
               "-o", os.devnull, # throw away the verbose output
               "--domtblout", "result.domtab",
               query_hmmfile]

    # Allow for disabling multithreading for HMMer3 calls in the command line
    if config.get('hmmer3') and 'multithreading' in config.hmmer3 and \
            not config.hmmer3.multithreading:
        command = command[0:1] + command[3:]

    with TemporaryDirectory(change=True):
        try:
            if use_tempfile:
                with open("input.fa", 'w') as handle:
                    handle.write(target_sequence)
                command.append("input.fa")
                run_result = execute(command)
            else:
                command.append('-')
                run_result = execute(command, stdin=target_sequence)
        except OSError:
            return []
        if not run_result.successful():
            logging.error('hmmsearch returned %d: %s while searching %s',
                          run_result.return_code, run_result.stderr, query_hmmfile)
            raise RuntimeError("Running hmmsearch failed.")
        return list(SearchIO.parse("result.domtab", 'hmmsearch3-domtab'))

def run_hmmpress(hmmfile):
    "Run hmmpress"
    command = ['hmmpress', hmmfile]
    try:
        run_result = execute(command)
        out = run_result.stdout
        err = run_result.stderr
        retcode = run_result.return_code
    except OSError as excep:
        retcode = 1
        err = str(excep)
        out = None
    return out, err, retcode

def run_hmmpfam2(query_hmmfile, target_sequence): # TODO cleanup
    "Run hmmpfam2"
    config = get_config()
    command = ["hmmpfam2", "--cpu", str(config.cpus),
               query_hmmfile, '-']

    # Allow to disable multithreading for HMMer2 calls in the command line #TODO fix options for this
    if config.get('hmmer2') and 'multithreading' in config.hmmer2 and \
            not config.hmmer2.multithreading:
        command = command[0:1] + command[3:]

    result = execute(command, stdin=target_sequence)
    if not result.successful():
        logging.debug('hmmpfam2 returned %d: %r while searching %r', result.return_code,
                        result.stderr, query_hmmfile)
        raise RuntimeError("hmmpfam2 problem while running %s", command)
    res_stream = StringIO(result.stdout)
    results = list(SearchIO.parse(res_stream, 'hmmer2-text'))
    return results

def run_fimo_simple(query_motif_file, target_sequence): # TODO cleanup
    "Run FIMO"
    command = ["fimo", "--text", "--verbosity", "1", query_motif_file, target_sequence]
    result = execute(command)
    if not result.successful():
        logging.debug('FIMO returned %d: %r while searching %r', result.return_code,
                        result.stderr, query_motif_file)
        raise RuntimeError("FIMO problem while running %s... %s", command, result.stderr[-100:])
    return result.stdout

def run_hmmscan(target_hmmfile, query_sequence, opts=None, results_file=None):
    "Run hmmscan on the inputs and return a list of QueryResults"
    if not query_sequence:
        raise ValueError("Cannot run hmmscan on empty sequence")

    config = get_config()
    command = ["hmmscan", "--cpu", str(config.cpus), "--nobias"]

    # Allow to disable multithreading for HMMer3 calls in the command line
    if config.get('hmmer3') and 'multithreading' in config.hmmer3 and \
            not config.hmmer3.multithreading: # TODO: ensure working
        command = command[0:1] + command[3:]

    if opts is not None:
        command.extend(opts)
    command.extend([target_hmmfile, '-'])
    result = execute(command, stdin=query_sequence)
    if not result.successful():
        raise RuntimeError('hmmscan returned %d: %r while scanning %r' % (
                           result.return_code, result.stderr[-100:].replace("\n", ""),
                           query_sequence[:100]))
    if results_file is not None:
        with open(results_file, 'w') as fh:
            fh.write(result.stdout)

    return list(SearchIO.parse(StringIO(result.stdout), 'hmmer3-text'))
