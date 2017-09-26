# License: GNU Affero General Public License v3 or later
# A copy of GNU AGPL v3 should have been included in this software package in LICENSE.txt.

import cProfile
from io import StringIO
import logging
import pstats
import os
from datetime import datetime

from Bio import SeqIO

from antismash.config import loader, update_config
from antismash.common import deprecated, serialiser
from antismash.common.module_results import ModuleResults
from antismash.modules import tta, genefinding, hmm_detection, clusterblast, \
                              dummy, lanthipeptides, smcogs
from antismash.outputs import html, svg

__version__ = "5.0.0alpha"

def get_all_modules():
    return get_detection_modules() + get_analysis_modules() + get_output_modules()

def get_detection_modules():
    return [hmm_detection, genefinding]

def get_analysis_modules():
    return [smcogs, tta, lanthipeptides, clusterblast, dummy]

def get_output_modules():
    return [html]

def setup_logging(logfile=None, verbose=False, debug=False):
    "Set up the logging output"

    def new_critical(*args): #TODO: temporary to make alpha issues more obvious
        msg = "\033[1;33m{}\033[0m".format(args[0])
        print(msg%args[1:])
    logging.critical = new_critical

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


def verify_options(options, modules):
    errors = []
    for module in modules:
        try:
            logging.debug("Checking options for %s", module.__name__)
            errors.extend(module.check_options(options))
        except ValueError as err:
            errors.append(str(err))
    if not errors:
        return True
    logging.error("Incompatible options detected:")
    for error in errors:
        print(error) # still commandline args, so don't use logging
    return False

def detect_signature_genes(seq_record, options):
    "Detect different secondary metabolite clusters based on HMM signatures"
    logging.info('Looking for secondary metabolite cluster signatures')
    hmm_detection.detect_signature_genes(seq_record, options)

def run_detection_stage(record, options, detection_modules):
    # strip any existing antismash results first
    deprecated.strip_record(record)

    detection_results = {}
    for module in detection_modules:
        if module.is_enabled(options):
            detection_results[module.NAME] = module.run_on_record(record, options)
    return detection_results

def regenerate_results_for_record(record, options, modules, previous_result):
    for module in modules:
        section = previous_result.pop(module.__name__, None)
        results = None
        if section:
            logging.debug("Regenerating results for %s", module.__name__)
            results = module.check_previous_results(section, record, options)
            if not results:
                logging.debug("Results could not be generated for %s", module.__name__)
            else:
                assert isinstance(results, ModuleResults)
                previous_result[module.__name__] = results
    return previous_result

def analyse_record(record, options, modules, previous_result):
    module_results = regenerate_results_for_record(record, options, modules, previous_result)

    # try to run the given modules over the record
    logging.info("Analysing record: %s", record.id)
    for module in modules:
        logging.debug("Checking if %s should be run", module.__name__)
        if not module.is_enabled(options):
            continue
        results = module_results.get(module.__name__)
        logging.info("Running %s", module.__name__)
        results = module.run_on_record(record, results, options)
        assert isinstance(results, ModuleResults)
        module_results[module.__name__] = results


def prepare_output_directory(name):
    if not name:
        return
    if os.path.exists(name):
        if not os.path.isdir(name):
            raise RuntimeError("Output directory %s exists and is not a directory" % name)
        logging.debug("Reusing output directory: %s", name)
    else:
        logging.debug("Creating output directory: %s", name)
        os.mkdir(name)


def write_profiling_results(profiler, target):
    stream = StringIO()
    sortby = 'tottime'
    stats = pstats.Stats(profiler, stream=stream).sort_stats(sortby)
    stats.dump_stats(target + ".bin")
    stats.print_stats(.25) # limit to the more meaningful first 25%
    stats.print_callers(.25)
    try:
        path_to_remove = os.path.dirname(os.path.realpath(__file__)) + os.path.sep
        open(target, "w").write(stream.getvalue().replace(path_to_remove, ""))
        logging.info("Profiling report written to %s", target)
    except IOError:
        #if can't save to file, print to terminal, but only the head
        logging.debug("Couldn't open file to store profiling output")
        stream.truncate(0)
        stats.print_stats(20) #first 20 lines only
        print(stream.getvalue())


def write_outputs(results, options):
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
    combined_filename = os.path.join(options.output_dir, results.input_file)
    combined_filename = os.path.splitext(combined_filename)[0] + ".gbk"
    logging.debug("Writing final genbank file to '%s'", combined_filename)
    SeqIO.write(bio_records, combined_filename, "genbank")


def annotate_records(results):
    for record, record_results in zip(results.records, results.results):
        logging.debug("Annotating record %s with results from: %s", record.id,
                      ", ".join([name.split()[0].split('.')[-1] for name in record_results]))
        for module, result in record_results.items():
            logging.debug(" Adding results from %s", module)
            assert isinstance(result, ModuleResults), type(result)
            result.add_to_record(record)

def read_data(sequence_file, options):
    if sequence_file:
        records = deprecated.parse_input_sequence(sequence_file, options)
        return serialiser.AntismashResults(sequence_file.rsplit(os.sep, 1)[-1],
                              records, [{}] * len(records), __version__)

    logging.debug("Attempting to reuse previous results in: %s", options.reuse_results)
    with open(options.reuse_results) as handle:
        contents = handle.read()
        if not contents:
            raise ValueError("No results contained in file: %s" % options.reuse_results)
    results = serialiser.AntismashResults.from_file(options.reuse_results)
    # hacky bypass to set output dir #TODO work out alternate method
    options.__dict__["output_dir"] = os.path.abspath(os.path.splitext(results.input_file)[0])
    return results

def check_prerequisites(modules):
    for module in modules:
        logging.debug("Checking prerequisites for %s", module.__name__)
        res = module.check_prereqs()
        if res:
            raise RuntimeError("Module failing prerequisite check: %s %s" %(
                            module.__name__, res))

def list_plugins(modules):
    print("Available plugins")
    max_name = 0
    for module in modules:
        max_name = max(max_name, len(module.NAME))
    format_string = " %-{}s:  %s".format(max_name)
    for module in modules:
        print(format_string % (module.NAME, module.SHORT_DESCRIPTION))

def run_antismash(sequence_file, options, detection_modules=None,
                  analysis_modules=None):
    logfile = options.logfile if 'logfile' in options else None
    setup_logging(logfile=logfile, verbose=options.verbose,
                  debug=options.debug)


    if detection_modules is None:
        detection_modules = get_detection_modules()
    if analysis_modules is None:
        analysis_modules = get_analysis_modules()

    modules = detection_modules + analysis_modules

    if options.list_plugins:
        list_plugins(modules)
        return 0

    options.all_enabled_modules = [module for module in modules if module.is_enabled(options)]
    options = update_config(options)
    loader.update_config_from_file() # TODO move earlier to run_antismash?

    check_prerequisites(modules)
    if options.check_prereqs_only:
        print("All prerequisites satisfied")
        return 0

    # start up profiling if relevant
    if options.profile:
        profiler = cProfile.Profile()
        profiler.enable()

    # ensure the provided options are valid
    if not verify_options(options, analysis_modules + detection_modules):
        return 1


    # check that at least one module will run
    if not any(module.is_enabled(options) for module in detection_modules + analysis_modules):
        raise ValueError("No detection or analysis modules enabled")

    start_time = datetime.now()

    prepare_output_directory(options.output_dir)

    results = read_data(sequence_file, options)

    results.records = deprecated.pre_process_sequences(results.records, options, genefinding)
    for seq_record, previous_result in zip(results.records, results.results):
        # skip if we're not interested in it
        if seq_record.skip:
            continue
        run_detection_stage(seq_record, options, detection_modules)
        analyse_record(seq_record, options, analysis_modules, previous_result)


    # Write results
    # TODO: include status logging, zipping, etc
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
        write_profiling_results(profiler,
                        os.path.join(options.output_dir, "profiling_results"))

    end_time = datetime.now()
    running_time = end_time - start_time

    logging.debug("antiSMASH calculation finished at %s; runtime: %s", str(end_time), str(running_time))

    # TODO: utils.log_status("antiSMASH status: SUCCESS")
    logging.debug("antiSMASH status: SUCCESS")
    return 0
