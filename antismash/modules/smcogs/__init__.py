# License: GNU Affero General Public License v3 or later
# A copy of GNU AGPL v3 should have been included in this software package in LICENSE.txt.

import logging
import os

from antismash.common import deprecated, path, subprocessing, hmmscan_refinement
from antismash.common.module_results import ModuleResults
from antismash.config.args import ModuleArgs

from .trees import generate_trees
from .classify import classify_genes, load_cog_annotations, write_smcogs_file

NAME = "smcogs"
SHORT_DESCRIPTION = NAME.capitalize()


def check_options(options):
    errors = []
    if options.smcogs_trees and (options.minimal and not options.smcogs_enabled):
        logging.debug("SMCOG trees enabled, but not classifications, running both anyway")
    return errors


def get_arguments():
    group = ModuleArgs("Basic analysis options", "smcogs", basic_help=True,
                       enabled_by_default=True)
    group.add_analysis_toggle('--smcogs-trees',
                              dest='smcogs_trees',
                              action='store_true',
                              default=False,
                              help="Generate phylogenetic trees of sec. "
                                   "met. cluster orthologous groups.")
    return group


def is_enabled(options):
    return not options.minimal or options.smcogs_enabled or options.smcogs_trees


def regenerate_previous_results(results, record, options):
    if not results or record.id != results["record_id"]:
        return None
    if options.smcogs_trees and not results["tree_paths"]:
        # trees have to be regenerated, so don't reuse
        logging.debug("Trees require recalculation")
        return None
    parsed = SMCOGResults.from_json(results)
    for tree_filename in parsed.tree_images.values():
        if not os.path.exists(os.path.join(parsed.relative_tree_path, tree_filename)):
            logging.debug("Tree image files missing and must be regenerated")
            return None
    return parsed


def check_prereqs():
    "Check if all required applications are around"
    failure_messages = []
    for binary_name in ['muscle', 'hmmscan', 'hmmpress', 'fasttree', 'java']:
        if path.locate_executable(binary_name) is None:
            failure_messages.append("Failed to locate file: %r" % binary_name)

    for hmm in ['smcogs.hmm']:
        hmm = path.get_full_path(__file__, 'data', hmm)
        if path.locate_file(hmm) is None:
            failure_messages.append("Failed to locate file %r" % hmm)
            continue
        for ext in ['.h3f', '.h3i', '.h3m', '.h3p']:
            binary = "%s%s" % (hmm, ext)
            if path.locate_file(binary) is None:
                # regenerate them
                result = subprocessing.run_hmmpress(hmm)
                if not result.successful():
                    failure_messages.append("Failed to hmmpress %s: %s" % (hmm, result.stderr.rstrip()))
                break
    return failure_messages


class SMCOGResults(ModuleResults):
    schema_version = 1

    def __init__(self, record_id):
        super().__init__(record_id)
        self.tree_images = {}  # gene_id -> tree filename
        self.best_hits = {}  # gene_id -> best HMM result
        self.relative_tree_path = None  # path where tree images are saved

    def to_json(self):
        return {"schema_version": self.schema_version,
                "record_id": self.record_id,
                "tree_paths": self.tree_images,
                "best_hits": {key: [getattr(val, attr) for attr in val.__slots__] for key, val in self.best_hits.items()},
                "image_dir": self.relative_tree_path}

    @staticmethod
    def from_json(json):
        if json.get("schema_version") != SMCOGResults.schema_version:
            return None
        results = SMCOGResults(json["record_id"])
        for hit, parts in json["best_hits"].items():
            results.best_hits[hit] = hmmscan_refinement.HMMResult(*parts)

        if json.get("image_dir"):
            results.relative_tree_path = json["image_dir"]
            results.tree_images = json["tree_paths"]
        return results

    def add_to_record(self, record):
        """ Annotate smCOGS in CDS features """
        functions = load_cog_annotations()
        logging.debug("annotating genes with SMCOGS info: %d genes", len(self.best_hits))
        for feature in record.get_cds_features_within_clusters():
            gene_id = feature.get_name()
            result = self.best_hits.get(gene_id)
            if result:  # TODO convert to qualifier like SecMetQualifier
                smcog_id, name = result.hit_id.split(':')
                feature.gene_functions.add(functions[smcog_id], "smcogs", "%s (Score: %g; E-value: %g)" % (result.hit_id, result.bitscore, result.evalue))
            if gene_id in self.tree_images:
                feature.notes.append("smCOG tree PNG image: smcogs/%s" % self.tree_images[gene_id])


def run_on_record(record, results, options):
    relative_output_dir = os.path.relpath(os.path.join(options.output_dir, "smcogs"), os.getcwd())
    smcogs_dir = os.path.abspath(relative_output_dir)
    if not os.path.exists(smcogs_dir):
        os.mkdir(smcogs_dir)

    if not results:
        results = SMCOGResults(record.id)

        genes = record.get_cds_features_within_clusters()
        hmm_results = classify_genes(genes)
        for gene in genes:
            gene_name = gene.get_name()
            hits = hmm_results.get(gene_name)
            if not hits:
                continue
            results.best_hits[gene.get_name()] = hits[0]
        write_smcogs_file(hmm_results, genes, deprecated.get_pksnrps_cds_features(record), options)

    if not results.tree_images and options.smcogs_trees:
        # create the smcogs output directory if required
        results.relative_tree_path = relative_output_dir
        original_dir = os.getcwd()
        os.chdir(smcogs_dir)  # TODO make a context manager
        nrpspks_genes = deprecated.get_pksnrps_cds_features(record)
        nrpspks_genes = []
        results.tree_images = generate_trees(smcogs_dir, hmm_results, genes, nrpspks_genes, options)

        os.chdir(original_dir)

    return results
