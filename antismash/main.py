# License: GNU Affero General Public License v3 or later
# A copy of GNU AGPL v3 should have been included in this software package in LICENSE.txt.

""" The core setup functions wrapping all detection, analysis, and output
    modules together.

    The intended entry point of antismash is run_antismash() in this file.
"""

from collections import defaultdict
import cProfile
from datetime import datetime
from io import StringIO
import logging
import os
import pstats
import shutil
import time
import tempfile
from types import ModuleType
from typing import Dict, List, Optional, Union

from Bio import SeqIO

from antismash.config import update_config
from antismash.common import serialiser, record_processing
from antismash.common.module_results import ModuleResults
from antismash.common.secmet import Record
from antismash.detection import genefinding, hmm_detection, nrps_pks_domains, full_hmmer, \
                                cassis, clusterfinder, cluster_hmmer
from antismash.modules import tta, clusterblast, lanthipeptides, smcogs, dummy, \
                              nrps_pks, thiopeptides, sactipeptides, lassopeptides, active_site_finder
from antismash.outputs import html, svg

__version__ = "5.0.0alpha"


def get_all_modules() -> List[ModuleType]:
    """ Return a list of default modules

        Arguments:
            None

        Returns:
            a list of modules
    """
    return get_detection_modules() + get_analysis_modules() + get_output_modules()


def get_detection_modules() -> List[ModuleType]:
    """ Return a list of default detection modules

        Arguments:
            None

        Returns:
            a list of modules
    """
    return [genefinding, hmm_detection, nrps_pks_domains, full_hmmer, cassis, clusterfinder,
            cluster_hmmer]


def get_analysis_modules() -> List[ModuleType]:
    """ Return a list of default analysis modules

        Arguments:
            None

        Returns:
            a list of modules
    """
    return [smcogs, tta, lanthipeptides, thiopeptides, nrps_pks, clusterblast,
            sactipeptides, lassopeptides, active_site_finder, dummy]


def get_output_modules() -> List[ModuleType]:
    """ Return a list of default output modules

        Arguments:
            None

        Returns:
            a list of modules
    """
    return [html]


def setup_logging(logfile: str = None, verbose: bool = False, debug: bool = False) -> None:
    """ Define the logging format, levels and outputs

        Arguments:
            logfile: None or the path to a file to write logging messages to
            verbose: whether to show INFO level messages and above
            debug: whether to show DEBUG level messages and above

        Returns:
            None
    """

    def new_critical(*args):  # TODO: temporary to make alpha issues more obvious
        """ make critical messages yellow and without the normal timestamp """
        msg = "\033[1;33m{}\033[0m".format(args[0])
        print(msg % args[1:])
    logging.critical = new_critical  # type: ignore

    log_level = logging.WARNING
    if debug:
        log_level = logging.DEBUG
    elif verbose:
        log_level = logging.INFO

    log_format = '%(levelname)-8s %(asctime)s   %(message)s'
    logging.basicConfig(format=log_format, level=log_level, datefmt="%d/%m %H:%M:%S")

    if not logfile:
        return

    if not (os.path.dirname(logfile) == "" or os.path.exists(os.path.dirname(logfile))):
        os.mkdir(os.path.dirname(logfile))
    handler = logging.FileHandler(logfile)
    if log_level == logging.WARNING:
        log_level = logging.INFO
    handler.setLevel(log_level)
    handler.setFormatter(logging.Formatter(fmt=log_format, datefmt="%d/%m %H:%M:%S"))
    logging.getLogger('').addHandler(handler)


def verify_options(options, modules: List[ModuleType]) -> bool:
    """ Find and display any incompatibilities in provided options

        Arguments:
            options: the options to check
            modules: the modules to check the options of

        Returns:
            True if no problems detected, otherwise False
    """
    errors = []  # type: List[str]
    for module in modules:
        try:
            logging.debug("Checking options for %s", module.__name__)
            errors.extend(module.check_options(options))
        except ValueError as err:
            errors.append(str(err))
    if not errors:
        return True

    logging.error("Incompatible options detected:\n  %s", "\n  ".join(errors))
    for error in errors:
        print(error)  # still commandline args, so don't use logging
    return False


def run_detection(record: Record, options, previous_result: Dict[str, Union[Dict, ModuleResults]]) -> Dict[str, float]:
    """ Detect different secondary metabolite clusters, PFAMs, and domains.

        Arguments:
            record: the Record to run detection over
            options: antiSMASH config
            previous_result: a dictionary mapping a module's name to results from
                             a previous run on this module, either as a
                             ModuleResults subclass or a JSON-like dictionary if
                             the module created it is no longer present

        Returns:
            the time taken by each detection module as a dictionary
    """
    # strip any existing antismash results first  # TODO: don't strip detection stage results if reusing
    record_processing.strip_record(record)

    timings = {}

    module_results = regenerate_results_for_record(record, options, get_detection_modules(),
                                                   previous_result)

    # run full genome detections
    for module in [full_hmmer]:
        run_module(record, module, options, module_results, timings)
        results = module_results.get(module.__name__)
        if results:
            results.add_to_record(record)

    # generate cluster predictions
    logging.info("Detecting secondary metabolite clusters")
    predictions = []
    for module in [hmm_detection, cassis, clusterfinder]:
        run_module(record, module, options, module_results, timings)
        results = module_results.get(module.__name__)
        if results:
            predictions.extend(results.get_predictions())

    # create merged clusters
    record.create_clusters_from_borders(predictions)
    for cluster in record.get_clusters():
        cluster.trim_overlapping()

    if not record.get_clusters():
        logging.debug("No clusters detected, skipping record")
        record.skip = "No clusters detected"
        return None

    logging.info("%d cluster(s) detected in record", len(record.get_clusters()))

    # finally, run any detection limited to genes in clusters
    for module in [nrps_pks_domains, cluster_hmmer]:
        run_module(record, module, options, module_results, timings)
        results = module_results.get(module.__name__)
        if results:
            results.add_to_record(record)

    return timings


def regenerate_results_for_record(record: Record, options, modules: List[ModuleType],
                                  previous_result: Dict[str, Dict]
                                  ) -> Dict[str, Optional[ModuleResults]]:
    """ Converts a record's JSON results to ModuleResults per module

        Arguments:
            record: the record to regenerate results for
            options: antismash Config
            modules: the modules to regenerate results of
            previous_result: a dict of the json results to convert, in the form:
                    {modulename : {module details}}

        Returns:
            the previous_result dict, with the values of all modules provided
            as an instance of ModuleResults or None if results don't apply or
            could not be regenerated
    """
    # skip if nothing to work with
    if not previous_result:
        return previous_result

    logging.debug("Regenerating results for record %s", record.id)
    for module in modules:
        section = previous_result.pop(module.__name__, None)
        results = None
        if section:
            logging.debug("Regenerating results for module %s", module.__name__)
            results = module.regenerate_previous_results(section, record, options)
            if not results:
                logging.debug("Results could not be generated for %s", module.__name__)
            else:
                assert isinstance(results, ModuleResults)
                previous_result[module.__name__] = results
    return previous_result


def run_module(record, module, options, module_results, timings) -> None:
    """ Run analysis modules on a record

        Arguments:
            record: the record to run the analysis on
            options: antismash Config
            modules: the modules to analyse with
                        each module will run only if enabled and not reusing all
                        results
            module_results: a dictionary of module name to json results,
                                json results will be replaced by ModuleResults
                                instances
            timings: a dictionary mapping module name to time taken for that
                     module, will be updated with the module timing

        Returns:
            the time taken to run the module as a float
        """

    logging.debug("Checking if %s should be run", module.__name__)
    if not module.is_enabled(options):
        return

    results = module_results.get(module.__name__)
    logging.info("Running %s", module.__name__)

    start = time.time()
    results = module.run_on_record(record, results, options)
    duration = time.time() - start

    assert isinstance(results, ModuleResults), "%s returned %s" % (module.__name__, type(results))
    module_results[module.__name__] = results
    timings[module.__name__] = duration


def analyse_record(record, options, modules, previous_result) -> Dict[str, float]:
    """ Run analysis modules on a record

        Arguments:
            record: the record to run the analysis on
            options: antismash Config
            modules: the modules to analyse with
                        each module will run only if enabled and not reusing all
                        results
            previous_result: a dictionary of module name to json results,
                                json results will be replaced by ModuleResults
                                instances

        Returns:
            a dictionary mapping module name to time taken
    """
    module_results = regenerate_results_for_record(record, options, modules, previous_result)
    timings = {}
    # try to run the given modules over the record
    logging.info("Analysing record: %s", record.id)
    for module in modules:
        run_module(record, module, options, module_results, timings)
    return timings


def prepare_output_directory(name: str, input_file: str) -> None:
    """ Ensure the ouptut directory exists and is usable

        Raises an exception if the directory is unusable

        Arguments:
            name: the path of the directory
            input_file: the path of the input file

        Returns:
            None
    """
    # if not supplied, set the output directory to be the sequence name
    if not name:
        name = os.path.abspath(os.path.splitext(os.path.basename(input_file))[0])
        update_config({"output_dir": name})

    if os.path.exists(name):
        if not os.path.isdir(name):
            raise RuntimeError("Output directory %s exists and is not a directory" % name)
        logging.debug("Reusing output directory: %s", name)
    else:
        logging.debug("Creating output directory: %s", name)
        os.mkdir(name)


def write_profiling_results(profiler, target) -> None:
    """ Write profiling files to file in human readable form and as a binary
        blob for external tool use (with the extra extension '.bin').

        If the file cannot be opened or written to, a shortened form will be
        written to stdout to avoid losing the data.

        Arguments:
            profiler: the profiler instance to log results of
            target: the path of the file to store reuslts in

        Returns:
            None
    """
    stream = StringIO()
    sortby = 'tottime'
    stats = pstats.Stats(profiler, stream=stream).sort_stats(sortby)
    stats.dump_stats(target + ".bin")
    stats.print_stats(.25)  # limit to the more meaningful first 25%
    stats.print_callers(.25)
    try:
        path_to_remove = os.path.dirname(os.path.realpath(__file__)) + os.path.sep
        open(target, "w").write(stream.getvalue().replace(path_to_remove, ""))
        logging.info("Profiling report written to %s", target)
    except IOError:
        # if can't save to file, print to terminal, but only the head
        logging.debug("Couldn't open file to store profiling output")
        stream.truncate(0)
        stats.print_stats(20)  # first 20 lines only
        print(stream.getvalue())


def write_outputs(results, options) -> None:
    """ Write output files (webpage, genbank files, etc) to the output directory

        Arguments:
            results: a serialiser.AntismashResults instance
            options: an antismash Config instance

        Returns:
            None
    """
    logging.debug("Creating results page")
    html.write(results.records, results.results, options)

    logging.debug("Creating results SVGs")
    svg.write(options, results.results)

    # convert records to biopython
    bio_records = [record.to_biopython() for record in results.records]

    logging.debug("Writing cluster-specific genbank files")
    for record, bio_record in zip(results.records, bio_records):
        for cluster in record.get_clusters():
            cluster.write_to_genbank(directory=options.output_dir, record=bio_record)

    # write records to an aggregate output
    base_filename = os.path.splitext(os.path.join(options.output_dir, results.input_file))[0]
    combined_filename = base_filename + ".gbk"
    logging.debug("Writing final genbank file to '%s'", combined_filename)
    SeqIO.write(bio_records, combined_filename, "genbank")

    zipfile = base_filename + ".zip"
    logging.debug("Zipping output to '%s'", zipfile)
    if os.path.exists(zipfile):
        os.remove(zipfile)

    with tempfile.NamedTemporaryFile(prefix="as_zip_tmp", suffix=".zip") as temp:
        shutil.make_archive(temp.name.replace(".zip", ""), "zip", root_dir=options.output_dir)
        shutil.copy(temp.name, zipfile)
    assert os.path.exists(zipfile)


def annotate_records(results) -> None:
    """ Annotates all analysed records with the results generated from them

        Arguments:
            results: a serialiser.AntismashResults instance

        Returns:
            None
    """
    detection_module_names = {mod.__name__ for mod in get_detection_modules()}
    for record, record_results in zip(results.records, results.results):
        if record.skip:
            logging.debug("Not annotating skipped record %s: %s", record.id, record.skip)
            continue
        if not record_results:
            logging.debug("No results for record %s, not annotating", record.id)
            continue
        logging.debug("Annotating record %s with results from: %s", record.id,
                      ", ".join([name.split()[0].split('.')[-1] for name in record_results]))
        for module, result in record_results.items():
            if module in detection_module_names:
                continue
            logging.debug(" Adding results from %s", module)
            assert isinstance(result, ModuleResults), type(result)
            result.add_to_record(record)


def read_data(sequence_file, options) -> serialiser.AntismashResults:
    """ Reads in the data to be used in the analysis run. Can be provided as
        as a sequence file (fasta/genbank) or as file of prior results

        Arguments:
            sequence_file: A fasta/genbank file to read (or None)
            options: An antismash Config instance

        Returns:
            a AntismashResults instance, populated only if reusing results

    """
    if not sequence_file and not options.reuse_results:
        raise ValueError("No sequence file or prior results to read")

    if sequence_file:
        records = record_processing.parse_input_sequence(sequence_file, options.taxon,
                                options.minlength, options.start, options.end)
        return serialiser.AntismashResults(sequence_file.rsplit(os.sep, 1)[-1],
                                           records, [{} for i in range(len(records))],
                                           __version__)

    logging.debug("Attempting to reuse previous results in: %s", options.reuse_results)
    with open(options.reuse_results) as handle:
        contents = handle.read()
        if not contents:
            raise ValueError("No results contained in file: %s" % options.reuse_results)
    results = serialiser.AntismashResults.from_file(options.reuse_results, options.taxon)
    return results


def check_prerequisites(modules: List[ModuleType]) -> None:
    """ Checks that each module's prerequisites are satisfied. If not satisfied,
        a RuntimeError is raised.

        Any issues found are logged to logging.error

        Arguments:
            modules: the modules to check

        Returns:
            None
    """
    errors_by_module = {}
    for module in modules:
        logging.debug("Checking prerequisites for %s", module.__name__)
        res = module.check_prereqs()
        if res:
            errors_by_module[module.__name__] = res
    if errors_by_module:
        for module, errors in errors_by_module.items():
            for error in errors:
                logging.error("%s: preqrequisite failure: %s", module, error)
        raise RuntimeError("Modules failing prerequisites")


def list_plugins(modules: List[ModuleType]) -> None:
    """ Prints the name and short description of the given modules

        Arguments:
            modules: the modules to display info of

        Returns:
            None
    """
    print("Available plugins")
    max_name = 0
    for module in modules:
        max_name = max(max_name, len(module.NAME))
    format_string = " %-{}s:  %s".format(max_name)
    for module in modules:
        print(format_string % (module.NAME, module.SHORT_DESCRIPTION))


def log_module_runtimes(timings: Dict[str, Dict[str, float]]) -> None:
    """ Log the aggregate time taken per module.

        Arguments:
            timings: a dictionary mapping record id to
                        a dictionary mapping module name to time taken

        Returns:
            None
    """
    total_times = defaultdict(lambda: 0.)  # type: Dict[str, float]
    for result in timings.values():
        for module, runtime in result.items():
            total_times[module] += runtime
    if not total_times:
        return
    logging.debug("Total times taken by modules")
    for module, runtime in sorted(total_times.items()):
        logging.debug("  %s: %.1fs", module, runtime)


def run_antismash(sequence_file: Optional[str], options) -> int:
    """ The complete antismash pipeline. Reads in data, runs detection and
        analysis modules over any records found, then outputs the results to
        file.

        Arguments:
            sequence_file: the sequence file to read in records from, can be
                            None if reusing results
            options: command line options as an argparse.Namespace
            detection_modules: None or a list of modules to use for detection,
                                if None defaults will be used
            analysis_modules: None or a list of modules to use for analysis,
                                if None defaults will be used

        Returns:
            0 if requested operations completed succesfully, otherwise 1
            Exceptions may also be raised
    """
    logfile = options.logfile
    setup_logging(logfile=logfile, verbose=options.verbose,
                  debug=options.debug)

    detection_modules = get_detection_modules()
    analysis_modules = get_analysis_modules()
    modules = detection_modules + analysis_modules

    if options.list_plugins:
        list_plugins(modules)
        return 0

    options.all_enabled_modules = [module for module in modules if module.is_enabled(options)]
    # converts from a namespace to an antismash.config.Config instance so
    # modules can't fiddle with it
    options = update_config(options)

    if options.check_prereqs_only:
        try:
            check_prerequisites(modules)
        except RuntimeError:
            print("Some module prerequisites not satisfied")
            return 1
        print("All prerequisites satisfied")
        return 0
    else:
        check_prerequisites(options.all_enabled_modules)

    # start up profiling if relevant
    if options.profile:
        profiler = cProfile.Profile()
        profiler.enable()

    # ensure the provided options are valid
    if not verify_options(options, options.all_enabled_modules):
        return 1  # TODO: change to a raise?

    # check that at least one module will run
    if not options.all_enabled_modules:
        raise ValueError("No detection or analysis modules enabled")

    start_time = datetime.now()

    results = read_data(sequence_file, options)

    # reset module timings
    results.timings_by_record = {}

    prepare_output_directory(options.output_dir, sequence_file or options.reuse_results)

    results.records = record_processing.pre_process_sequences(results.records,
                                                              options, genefinding)
    for seq_record, previous_result in zip(results.records, results.results):
        # skip if we're not interested in it
        if seq_record.skip:
            continue
        timings = run_detection(seq_record, options, previous_result)
        # and skip analysis if detection didn't find anything
        if not seq_record.get_clusters():
            continue
        analysis_timings = analyse_record(seq_record, options, analysis_modules, previous_result)
        timings.update(analysis_timings)
        results.timings_by_record[seq_record.id] = timings

    # Write results
    json_filename = os.path.join(options.output_dir, results.input_file)
    json_filename = os.path.splitext(json_filename)[0] + ".json"
    logging.debug("Writing json results to '%s'", json_filename)
    results.write_to_file(json_filename)

    # now that the json is out of the way, annotate the record
    # otherwise we could double annotate some areas
    annotate_records(results)

    # create relevant output files
    write_outputs(results, options)

    # save profiling data
    if options.profile:
        profiler.disable()
        write_profiling_results(profiler, os.path.join(options.output_dir,
                                                       "profiling_results"))

    running_time = datetime.now() - start_time

    # display module runtimes before total time
    if options.debug:
        log_module_runtimes(results.timings_by_record)

    logging.debug("antiSMASH calculation finished at %s; runtime: %s",
                  str(datetime.now()), str(running_time))

    logging.info("antiSMASH status: SUCCESS")
    return 0
