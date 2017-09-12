# License: GNU Affero General Public License v3 or later
# A copy of GNU AGPL v3 should have been included in this software package in LICENSE.txt.

import json
import logging
import os
from Bio import SeqIO
from datetime import datetime

from antismash.config import loader
from antismash.common import deprecated, serialiser
from antismash.common.module_results import ModuleResults
from antismash.common.secmet import Record
from antismash.modules import tta, genefinding, hmm_detection, clusterblast, \
                              dummy, lanthipeptides, smcogs
from antismash.outputs import html, svg

def get_all_modules():
    return get_detection_modules() + get_analysis_modules() + get_output_modules()

def get_detection_modules():
    return [hmm_detection, genefinding]

def get_analysis_modules():
    return [smcogs, tta, clusterblast, lanthipeptides, dummy]

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
    fh = logging.FileHandler(logfile)
    if log_level == logging.WARNING:
        log_level = logging.INFO
    fh.setLevel(log_level)
    fh.setFormatter(logging.Formatter(fmt=log_format, datefmt="%d/%m %H:%M:%S"))
    logging.getLogger('').addHandler(fh)


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

def analyse_record(record, options, modules, previous_result):
    # ensure a minimum level of result format
    if not previous_result:
        previous_result["record_id"] = record.id
        previous_result["modules"] = {}

    # check that at leats one module will run
    if not any(module.is_enabled(options) for module in modules):
        logging.info("Skipping record, no modules enabled for: %s", record.id)
        return False

    # try to run the given modules over the record
    logging.info("Analysing record: %s", record.id)
    for module in modules:
        logging.debug("Checking if %s should be run", module.__name__)
        section = previous_result.get("modules", {}).get(module.__name__)
        results = module.check_previous_results(section, record, options)
        assert results is None or isinstance(results, ModuleResults)
        if results:
            if module.is_enabled(options):
                logging.info("Skipping %s, reusing previous results", module.__name__)
        elif module.is_enabled(options):
            logging.info("Running %s", module.__name__)
            # we've checked before, but make sure the options make sense
            assert not module.check_options(options) # TODO: change to return truthy if good
            results = module.run_on_record(record, options)
            if "hmm_detection" not in module.__name__: # TODO: work out if we can keep these
                assert isinstance(results, ModuleResults)
        previous_result["modules"][module.__name__] = results


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


def run_antismash(sequence_file, options, detection_modules=None,
                  analysis_modules=None):
    setup_logging(logfile=options.get('logfile', None), verbose=options.verbose,
                  debug=options.debug)
    loader.update_config_from_file()

    if detection_modules is None:
        detection_modules = get_detection_modules()
    if analysis_modules is None:
        analysis_modules = get_analysis_modules()

    # ensure the provided options are valid
    if not verify_options(options, analysis_modules + detection_modules):
        return 1

    for module in detection_modules + analysis_modules:
        logging.debug("Checking prerequisites for %s", module.__name__)
        res = module.check_prereqs()
        if res:
            raise RuntimeError("Module failing prerequisite check: %s %s" %(
                            module.__name__, res))

    start_time = datetime.now()

    prepare_output_directory(options.output_dir)

    if sequence_file:
        seq_records = deprecated.parse_input_sequence(sequence_file, options)
        results = [{}] * len(seq_records)
    else:
        logging.debug("Attempting to reuse previous results in: %s", options.reuse_results)
        with open(options.reuse_results) as handle:
            contents = handle.read()
            if not contents:
                raise ValueError("No results contained in file: %s" % options.reuse_results)
        results = json.loads(contents) #TODO clean up
        seq_records = serialiser.read_records(options.reuse_results)
        seq_records = [Record.from_biopython(record) for record in seq_records]

    seq_records = deprecated.pre_process_sequences(seq_records, options, genefinding)
    for seq_record, previous_result in zip(seq_records, results):
        # skip if we're not interested in it
        if seq_record.skip:
            continue
        run_detection_stage(seq_record, options, detection_modules)
        analyse_record(seq_record, options, analysis_modules, previous_result)


    # Write results
    logging.debug("Creating results page")
    html.write(seq_records, results, options)
    logging.debug("Creating results SVGs")
    svg.write(seq_records, options, results)
    # TODO: include status logging, zipping, etc
    logging.debug("Writing cluster-specific genbank files")
    for record in seq_records:
        for cluster in record.get_clusters():
            cluster.write_to_genbank(directory=options.output_dir)
    logging.debug("Writing genbank file to 'temp.gbk'")
    seq_records = [record.to_biopython() for record in seq_records]
    SeqIO.write(seq_records, "temp.gbk", "genbank")
    logging.debug("Writing json results to 'temp.json'")
    serialiser.write_records(seq_records, results, "temp.json")


    end_time = datetime.now()
    running_time = end_time - start_time

    logging.debug("antiSMASH calculation finished at %s; runtime: %s", str(end_time), str(running_time))

    # TODO: utils.log_status("antiSMASH status: SUCCESS")
    logging.debug("antiSMASH status: SUCCESS")
    return 0
