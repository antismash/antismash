# License: GNU Affero General Public License v3 or later
# A copy of GNU AGPL v3 should have been included in this software package in LICENSE.txt.

""" A collection of functions and related classes for executing external
    commands from within antismash.
"""

from io import StringIO
import logging
import multiprocessing
import sys
import os
from subprocess import PIPE, Popen, TimeoutExpired
from tempfile import NamedTemporaryFile
from typing import Any, Callable, Dict, IO, Iterable, List, Optional, Union
import warnings

from helperlibs.wrappers.io import TemporaryDirectory

from antismash.config import get_config

from .fasta import write_fasta, read_fasta

# Don't display the SearchIO experimental warning, we know this.
with warnings.catch_warnings():
    warnings.simplefilter("ignore")
    from Bio import SearchIO


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


def run_hmmsearch(query_hmmfile: str, target_sequence: str, use_tempfile: bool = False
                  ) -> List[SearchIO._model.query.QueryResult]:  # pylint: disable=protected-access
    """ Run hmmsearch on a HMM file and a fasta input

        Arguments:
            query_hmmfile: the path to the HMM file
            target_sequence: the fasta input to search as a string
            use_tempfile: if True, a tempfile will be written for the fasta input
                          instead of piping

        Returns:
            a list of hmmsearch results as parsed by SearchIO
    """
    config = get_config()
    command = ["hmmsearch", "--cpu", str(config.cpus),
               "-o", os.devnull,  # throw away the verbose output
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


def run_hmmpress(hmmfile: str) -> RunResult:
    """ Run hmmpress on a HMMer model, overwriting any previous generated files
        (e.g. '.h3i').

        Arguments:
            hmmfile: the path to the HMMer model

        Returns:
            a RunResult instance
    """
    return execute(["hmmpress", "-f", hmmfile])


def run_hmmpfam2(query_hmmfile: str, target_sequence: str, extra_args: List[str] = None
                 ) -> List[SearchIO._model.query.QueryResult]:  # pylint: disable=protected-access
    """ Run hmmpfam2 over the provided HMM file and fasta input

        Arguments:
            query_hmmfile: the HMM file to use
            target_sequence: a string in fasta format of the sequence to run

        Returns:
            a list of results as parsed by SearchIO
    """
    config = get_config()
    command = ["hmmpfam2"]

    # Allow to disable multithreading for HMMer2 calls in the command line #TODO fix options for this
    if config.get('hmmer2') and 'multithreading' in config.hmmer2 and \
            config.hmmer2.multithreading:
        command.extend(["--cpu", str(config.cpus)])
    if extra_args:
        command.extend(extra_args)
    command.extend([query_hmmfile, '-'])

    result = execute(command, stdin=target_sequence)
    if not result.successful():
        logging.debug('hmmpfam2 returned %d: %r while searching %r', result.return_code,
                      result.stderr, query_hmmfile)
        raise RuntimeError("hmmpfam2 problem while running %s: %s" % (command, result.stderr))
    res_stream = StringIO(result.stdout)
    return list(SearchIO.parse(res_stream, 'hmmer2-text'))


def run_fimo_simple(query_motif_file: str, target_sequence: str) -> str:
    """ Runs FIMO on the provided inputs

        Arguments:
            query_motif_file: the path to the file containing query motifs
            target_sequence: the path to the file containing input sequences

        Returns:
            the output from running FIMO
    """
    command = ["fimo", "--text", "--verbosity", "1", query_motif_file, target_sequence]
    result = execute(command)
    if not result.successful():
        logging.debug('FIMO returned %d: %r while searching %r', result.return_code,
                      result.stderr, query_motif_file)
        raise RuntimeError("FIMO problem while running %s... %s" % (command, result.stderr[-100:]))
    return result.stdout


def run_hmmscan(target_hmmfile: str, query_sequence: str, opts: List[str] = None,
                results_file: str = None) -> List[SearchIO._model.query.QueryResult]:  # pylint: disable=protected-access
    """ Runs hmmscan on the inputs and return a list of QueryResults

        Arguments:
            target_hmmfile: the path to a HMM file to use in scanning
            query_sequence: a string containing input sequences in fasta format
            opts: a list of extra arguments to pass to hmmscan, or None
            results_file: a path to keep a copy of hmmscan results in, if provided

        Returns:
            a list of QueryResults as parsed from hmmscan output by SearchIO

    """
    if not query_sequence:
        raise ValueError("Cannot run hmmscan on empty sequence")

    config = get_config()
    command = ["hmmscan", "--cpu", str(config.cpus), "--nobias"]

    # Allow to disable multithreading for HMMer3 calls in the command line
    if config.get('hmmer3') and 'multithreading' in config.hmmer3 and \
            not config.hmmer3.multithreading:  # TODO: ensure working
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


def run_muscle_single(seq_name: str, seq: str, comparison_file: str) -> Dict[str, str]:
    """ Runs muscle over a single sequence against a comparison file in profile
        mode and returns a dictionary of the resulting alignments

        Arguments:
            seq_name: the name of the query
            seq: the sequence to align
            comparison_file: the path of the file containing comparison sequences

        Returns:
            a dictionary mapping sequence name (query or reference) to alignment
    """
    with NamedTemporaryFile(mode="w+") as temp_in:
        with NamedTemporaryFile(mode="w+") as temp_out:
            write_fasta([seq_name], [seq], temp_in.name)
            # Run muscle and collect sequence positions from file
            result = execute(["muscle", "-profile", "-quiet",
                              "-in1", comparison_file,
                              "-in2", temp_in.name,
                              "-out", temp_out.name])
            if not result.successful():
                raise RuntimeError("muscle returned %d: %r while comparing query named %s" % (
                                   result.return_code, result.stderr.replace("\n", ""),
                                   seq_name))
            fasta = read_fasta(temp_out.name)
    return fasta


def run_blastp(target_blastp_database: str, query_sequence: str,
               opts: List[str] = None, results_file: str = None
               ) -> List[SearchIO._model.query.QueryResult]:
    """ Runs blastp over a single sequence against a database and returns the
        results as parsed by Bio.SearchIO.

        Arguments:
            target_blastp_database: the blastp database to compare to
            query_sequence: the sequence being compared
            opts: a list of extra arguments to pass to blastp, or None
            results_file: a path to keep a copy of blastp results in, if provided

        Returns:
            a list of QueryResults as parsed from blast output by SearchIO
    """
    if not query_sequence:
        raise ValueError("Cannot run blastp on empty sequence")

    config = get_config()
    command = ["blastp", "-num_threads", str(config.cpus), "-db", target_blastp_database]

    if opts is not None:
        command.extend(opts)

    result = execute(command, stdin=query_sequence)
    if not result.successful():
        raise RuntimeError('blastp returned %d: %r while scanning %r' % (
                           result.return_code, result.stderr.replace("\n", ""),
                           query_sequence[:100]))

    if results_file is not None:
        with open(results_file, 'w') as fh:
            fh.write(result.stdout)

    return list(SearchIO.parse(StringIO(result.stdout), 'blast-text'))


def run_diamond(query_file: str, database_file: str, mode: str = "blastp",
                opts: Optional[List[str]] = None) -> str:
    """ Runs diamond, comparing the given query to the given database

        Arguments:
            query_file: the path of query sequence file
            database_file: the path of the database to compare to
            mode: the mode to use (defaults to blastp)
            opts: any extra options to pass to diamond

        Returns:
            the output from running diamond
    """
    with TemporaryDirectory() as temp_dir:
        command = [
            "diamond",
            mode,
            "--db", database_file,
            "--threads", str(get_config().cpus),
            "--query", query_file,
            "--tmpdir", temp_dir,
        ]
        if opts:
            command.extend(opts)
        result = execute(command)
        if not result.successful():
            raise RuntimeError("diamond failed to run: %s -> %s" % (command, result.stderr[-100:]))
    return result.stdout
