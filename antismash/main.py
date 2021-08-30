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
import glob
import logging
import os
import pstats
import shutil
import time
import tempfile
from typing import cast, Any, Dict, List, Optional, Tuple, Union

from Bio import SeqIO
from Bio.SeqRecord import SeqRecord

from antismash.config import (
    ConfigType,
    get_config,
    update_config,
)
from antismash.common import logs, record_processing, serialiser
from antismash.common.module_results import ModuleResults, DetectionResults
from antismash.common.secmet import Record
from antismash.common import subprocessing
from antismash.detection import (cassis,
                                 cluster_hmmer,
                                 full_hmmer,
                                 genefinding,
                                 hmm_detection,
                                 nrps_pks_domains,
                                 genefunctions,
                                 sideloader,
                                 tigrfam,
                                 )
from antismash.modules import (active_site_finder,
                               clusterblast,
                               cluster_compare,
                               lanthipeptides,
                               lassopeptides,
                               nrps_pks,
                               pfam2go,
                               rrefinder,
                               sactipeptides,
                               smcog_trees,
                               t2pks,
                               thiopeptides,
                               tta,
                               )
from antismash.outputs import html, svg
from antismash.custom_typing import AntismashModule

__version__ = "6.0.2"


def get_all_modules() -> List[AntismashModule]:
    """ Return a list of default modules

        Arguments:
            None

        Returns:
            a list of modules
    """
    return get_detection_modules() + get_analysis_modules() + get_output_modules()


def get_detection_modules() -> List[AntismashModule]:
    """ Return a list of default detection modules

        Arguments:
            None

        Returns:
            a list of modules
    """
    return [genefinding, hmm_detection, nrps_pks_domains, full_hmmer, cassis,  # type: ignore
            cluster_hmmer, genefunctions, sideloader, tigrfam]  # type: ignore


def get_analysis_modules() -> List[AntismashModule]:
    """ Return a list of default analysis modules

        Arguments:
            None

        Returns:
            a list of modules
    """
    return [smcog_trees, tta, lanthipeptides, thiopeptides, nrps_pks, clusterblast,  # type: ignore
            sactipeptides, lassopeptides, active_site_finder, pfam2go, t2pks, # type: ignore
            rrefinder, cluster_compare,  # type: ignore
           ]


def get_output_modules() -> List[AntismashModule]:
    """ Return a list of default output modules

        Arguments:
            None

        Returns:
            a list of modules
    """
    return [html]  # type: ignore  # a lot of casting avoided


def verify_options(options: ConfigType, modules: List[AntismashModule]) -> bool:
    """ Find and display any incompatibilities in provided options

        Arguments:
            options: the options to check
            modules: the modules to check the options of

        Returns:
            True if no problems detected, otherwise False
    """
    errors: List[str] = []
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


def run_detection(record: Record, options: ConfigType,
                  module_results: Dict[str, Union[ModuleResults, Dict[str, Any]]]) -> Dict[str, float]:
    """ Detect different secondary metabolite clusters, PFAMs, and domains.

        Arguments:
            record: the Record to run detection over
            options: antiSMASH config
            module_results: a dictionary mapping a module's name to results from
                            a previous run on this module, as a ModuleResults subclass
                            or in JSON form

        Returns:
            the time taken by each detection module as a dictionary
    """
    timings: Dict[str, float] = {}

    # run full genome detections
    for module in [full_hmmer]:
        run_module(record, cast(AntismashModule, module), options, module_results, timings)
        results = module_results.get(module.__name__)
        if results:
            assert isinstance(results, ModuleResults)
            logging.debug("Adding detection results from %s to record", module.__name__)
            results.add_to_record(record)

    # generate cluster predictions
    logging.info("Detecting secondary metabolite clusters")
    for module in [sideloader, hmm_detection, cassis]:
        run_module(record, cast(AntismashModule, module), options, module_results, timings)
        results = module_results.get(module.__name__)
        if results:
            assert isinstance(results, DetectionResults)
            for protocluster in results.get_predicted_protoclusters():
                record.add_protocluster(protocluster)
            for region in results.get_predicted_subregions():
                record.add_subregion(region)

    logging.debug("%d protoclusters found", len(record.get_protoclusters()))
    logging.debug("%d subregions found", len(record.get_subregions()))

    record.create_candidate_clusters()
    record.create_regions()

    if not record.get_regions():
        logging.info("No regions detected, skipping record")
        record.skip = "No regions detected"
        return timings

    logging.info("%d region(s) detected in record", len(record.get_regions()))

    # finally, run any detection limited to genes in clusters
    for module in [nrps_pks_domains, cluster_hmmer, genefunctions, tigrfam]:
        run_module(record, cast(AntismashModule, module), options, module_results, timings)
        results = module_results.get(module.__name__)
        if results:
            assert isinstance(results, ModuleResults)
            logging.debug("Adding detection results from %s to record", module.__name__)
            results.add_to_record(record)

    return timings


def run_module(record: Record, module: AntismashModule, options: ConfigType,
               module_results: Dict[str, Union[ModuleResults, Dict[str, Any]]],
               timings: Dict[str, float]
               ) -> None:
    """ Run a module on a record

        Arguments:
            record: the record to run the analysis on
            module: the module to run, only run if enabled and not reusing results
            options: antismash Config
            module_results: a dictionary of module name to ModuleResults
                            instances or their JSON representations,
                            updated if the module runs
            timings: a dictionary mapping module name to time taken for that
                     module, will be updated with the module timing

        Returns:
            None
    """
    previous_results = module_results.pop(module.__name__, None)
    results = None
    if previous_results is not None:
        assert isinstance(previous_results, dict)
        logging.debug("Regenerating results for %s", module.__name__)
        results = module.regenerate_previous_results(previous_results, record, options)
        if results:
            module_results[module.__name__] = results
    assert results is None or isinstance(results, ModuleResults)

    logging.debug("Checking if %s should be run", module.__name__)
    if not module.is_enabled(options):
        return

    logging.info("Running %s", module.__name__)

    start = time.time()
    results = module.run_on_record(record, results, options)
    duration = time.time() - start

    assert isinstance(results, ModuleResults), "%s returned %s" % (module.__name__, type(results))
    module_results[module.__name__] = results
    timings[module.__name__] = duration


def analyse_record(record: Record, options: ConfigType, modules: List[AntismashModule],
                   previous_result: Dict[str, Union[ModuleResults, Dict[str, Any]]]) -> Dict[str, float]:
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
    timings: Dict[str, float] = {}
    # try to run the given modules over the record
    for module in modules:
        run_module(record, module, options, previous_result, timings)
    return timings


def prepare_output_directory(name: str, input_file: str) -> None:
    """ Ensure the ouptut directory exists and is usable

        Raises an exception if the directory is unusable,
        or if results not being reused and directory not empty

        Arguments:
            name: the path of the directory
            input_file: the path of the input file

        Returns:
            None
    """
    # if not supplied, set the output directory to be the sequence name
    input_prefix = os.path.basename(canonical_base_filename(input_file, "", get_config()))
    if not name:
        name = os.path.abspath(input_prefix)
        update_config({"output_dir": name})

    if os.path.exists(name):
        if not os.path.isdir(name):
            raise RuntimeError("Output directory %s exists and is not a directory" % name)
        # not empty (apart from a possible input dir), and not reusing its results
        if not input_file.endswith(".json") and \
                list(filter(_ignore_patterns, glob.glob(os.path.join(name, "*")))):
            raise RuntimeError("Output directory contains other files, aborting for safety")

        # --reuse
        logging.debug("Removing existing region genbank files")
        for genbank in glob.glob(os.path.join(name, "*.region???.gbk")):
            os.remove(genbank)
        logging.debug("Reusing output directory: %s", name)
    else:
        logging.debug("Creating output directory: %s", name)
        os.mkdir(name)


def _ignore_patterns(entry: str) -> bool:
    """File name patterns that we want to ignore for the "outdir is empty" check."""
    config = get_config()
    if entry.endswith('/input') and os.path.isdir(entry):
        return False
    if os.path.abspath(entry) == os.path.abspath(config.logfile):
        return False

    return True


def write_profiling_results(profiler: cProfile.Profile, target: str) -> None:
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


def add_antismash_comments(records: List[Tuple[Record, SeqRecord]], options: ConfigType) -> None:
    """ Add antismash meta-annotation to records for genbank output

        Arguments:
            records: a list of Record, SeqRecord pairs
            options: antismash options

        Returns:
            None

    """
    if not records:
        return
    shared = []
    # include start/end details if relevant
    if options.start != -1 or options.end != -1:
        start = 1 if options.start == -1 else options.start
        # start/end is only valid for single records, as per record_processing
        assert len(records) == 1
        end = len(records[0][0].seq) if options.end == -1 else options.end
        shared.append((
            "NOTE: This is an extract from the original record!\n"
            "Starting at  :: {start}\n"
            "Ending at    :: {end}\n"
        ).format(start=start, end=end))
    antismash_comment = (
        "##antiSMASH-Data-START##\n"
        "Version      :: {version}\n"
        "Run date     :: {date}\n"
        "%s"
        "##antiSMASH-Data-END##"
        ).format(
            version=options.version,
            date=str(datetime.now().strftime("%Y-%m-%d %H:%M:%S")),
        )
    for record, bio_record in records:
        extras: List[str] = []
        if record.original_id:
            extras.append("Original ID  :: %s\n" % record.original_id)

        comment = antismash_comment % "".join(extras + shared)
        if 'comment' in bio_record.annotations:
            bio_record.annotations['comment'] += '\n' + comment
        else:
            bio_record.annotations['comment'] = comment


def write_outputs(results: serialiser.AntismashResults, options: ConfigType) -> None:
    """ Write output files (webpage, genbank files, etc) to the output directory

        Arguments:
            results: a serialiser.AntismashResults instance
            options: an antismash Config instance

        Returns:
            None
    """
    # don't use results for which the module no longer exists to regenerate/calculate
    module_results_per_record = []
    for record_results in results.results:
        record_result = {}
        for module_name, result in record_results.items():
            if isinstance(result, ModuleResults):
                record_result[module_name] = result
        module_results_per_record.append(record_result)

    logging.debug("Creating results page")
    html.write(results.records, module_results_per_record, options)

    logging.debug("Creating results SVGs")
    svg.write(options, module_results_per_record)

    # convert records to biopython
    bio_records = [record.to_biopython() for record in results.records]

    # add antismash meta-annotation to records
    add_antismash_comments(list(zip(results.records, bio_records)), options)

    logging.debug("Writing cluster-specific genbank files")
    for record, bio_record in zip(results.records, bio_records):
        for region in record.get_regions():
            region.write_to_genbank(directory=options.output_dir, record=bio_record)

    # write records to an aggregate output
    base_filename = canonical_base_filename(results.input_file, options.output_dir, options)
    combined_filename = base_filename + ".gbk"
    logging.debug("Writing final genbank file to '%s'", combined_filename)
    SeqIO.write(bio_records, combined_filename, "genbank")

    zipfile = base_filename + ".zip"
    if os.path.exists(zipfile):
        os.remove(zipfile)
    if not options.skip_zip_file:
        logging.debug("Zipping output to '%s'", zipfile)
        with tempfile.NamedTemporaryFile(prefix="as_zip_tmp", suffix=".zip") as temp:
            shutil.make_archive(temp.name.replace(".zip", ""), "zip", root_dir=options.output_dir)
            shutil.copy(temp.name, zipfile)
        assert os.path.exists(zipfile)


def canonical_base_filename(input_file: str, directory: str, options: ConfigType) -> str:
    """Generate a canonical base filename if one isn't specified in the options."""
    if options.output_basename:
        base_filename = options.output_basename
    else:
        base_filename, ext = os.path.splitext(os.path.basename(input_file))
        if ext.lower() in (".gz", ".bz", ".xz"):
            base_filename, _ = os.path.splitext(base_filename)
        update_config({"output_basename": base_filename})

    return os.path.join(directory, base_filename)


def annotate_records(results: serialiser.AntismashResults) -> None:
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


def read_data(sequence_file: Optional[str], options: ConfigType) -> serialiser.AntismashResults:
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
                                options.minlength, options.start, options.end,
                                gff_file=options.genefinding_gff3)
        results = serialiser.AntismashResults(sequence_file.rsplit(os.sep, 1)[-1],
                                              records, [{} for i in records],
                                              __version__, taxon=options.taxon)
        update_config({"input_file": os.path.splitext(results.input_file)[1]})
    else:
        logging.debug("Attempting to reuse previous results in: %s", options.reuse_results)
        with open(options.reuse_results) as handle:
            contents = handle.read()
            if not contents:
                raise ValueError("No results contained in file: %s" % options.reuse_results)
        results = serialiser.AntismashResults.from_file(options.reuse_results)
        for record in results.records:
            record.strip_antismash_annotations()
        if options.taxon != results.taxon:
            logging.info("Reusing taxon %s from prior results", results.taxon)
            update_config({"taxon": results.taxon})

    update_config({"input_file": os.path.splitext(results.input_file)[0]})
    return results


def prepare_module_data(modules: Optional[List[AntismashModule]] = None) -> None:
    """ Calls all given modules' data preparation methods.

        Arguments:
            modules: a list of modules to use, if None all module will be used

        Returns:
            None
    """
    if modules is None:
        modules = get_all_modules()
    for module in modules:
        if hasattr(module, "prepare_data"):
            module.prepare_data()


def check_prerequisites(modules: List[AntismashModule], options: ConfigType) -> None:
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
        res = module.check_prereqs(options)
        if res:
            errors_by_module[module.__name__] = res
    if errors_by_module:
        for module_name, errors in errors_by_module.items():
            for error in errors:
                logging.error("%s: preqrequisite failure: %s", module_name, error)
        raise RuntimeError("Modules failing prerequisites")


def list_plugins(modules: List[AntismashModule]) -> None:
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
    total_times:Dict[str, float] = defaultdict(lambda: 0.)
    for result in timings.values():
        for module, runtime in result.items():
            total_times[module] += runtime
    if not total_times:
        return
    logging.debug("Total times taken by modules")
    for module, runtime in sorted(total_times.items()):
        logging.debug("  %s: %.1fs", module, runtime)


def run_antismash(sequence_file: Optional[str], options: ConfigType) -> int:
    """ The complete antismash pipeline. Reads in data, runs detection and
        analysis modules over any records found, then outputs the results to
        file.

        Arguments:
            sequence_file: the sequence file to read in records from, can be
                            None if reusing results
            options: command line options

        Returns:
            0 if requested operations completed succesfully, otherwise 1
            Exceptions may also be raised
    """

    with logs.changed_logging(logfile=options.logfile, verbose=options.verbose,
                              debug=options.debug):
        result = _run_antismash(sequence_file, options)
    return result


def _run_antismash(sequence_file: Optional[str], options: ConfigType) -> int:
    """ The real run_antismash, assumes logging is set up around it """
    logging.info("antiSMASH version: %s", options.version)
    _log_found_executables(options)

    detection_modules = get_detection_modules()
    analysis_modules = get_analysis_modules()
    output_modules = get_output_modules()
    modules = detection_modules + analysis_modules + output_modules

    if options.list_plugins:
        list_plugins(modules)
        return 0

    options.all_enabled_modules = list(filter(lambda x: x.is_enabled(options), modules))

    if options.check_prereqs_only:
        try:
            check_prerequisites(modules, options)
        except RuntimeError:
            print("Some module prerequisites not satisfied")
            return 1
        print("All prerequisites satisfied")
        return 0

    check_prerequisites(options.all_enabled_modules, options)

    # start up profiling if relevant
    if options.profile:
        profiler = cProfile.Profile()
        profiler.enable()

    # ensure the provided options are valid
    if not verify_options(options, options.all_enabled_modules):
        return 1

    # check that at least one module will run
    if not options.all_enabled_modules:
        raise ValueError("No detection or analysis modules enabled")

    start_time = datetime.now()

    results = read_data(sequence_file, options)

    # reset module timings
    results.timings_by_record.clear()

    prepare_output_directory(options.output_dir, sequence_file or options.reuse_results)

    results.records = record_processing.pre_process_sequences(results.records, options,
                                                              cast(AntismashModule, genefinding))
    for record, module_results in zip(results.records, results.results):
        # skip if we're not interested in it
        if record.skip:
            continue
        logging.info("Analysing record: %s", record.id)
        timings = run_detection(record, options, module_results)
        # and skip analysis if detection didn't find anything
        if not record.get_regions():
            continue
        analysis_timings = analyse_record(record, options, analysis_modules, module_results)
        timings.update(analysis_timings)
        results.timings_by_record[record.id] = timings

    # Write results
    json_filename = canonical_base_filename(results.input_file, options.output_dir, options)
    json_filename += ".json"
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
                  datetime.now().strftime("%Y-%m-%d %H:%M:%S"), str(running_time))

    logging.info("antiSMASH status: SUCCESS")
    return 0


def _log_found_executables(options: ConfigType) -> None:
    for binary, path in vars(options.executables).items():
        version = ""
        version_getter = getattr(subprocessing, "run_{}_version".format(binary), None)
        if callable(version_getter):
            # pylint doesn't seem to understand this
            version = " ({})".format(version_getter())  # pylint: disable=not-callable
        logging.info("%s using executable: %s%s", binary, path, version)
