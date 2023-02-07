# License: GNU Affero General Public License v3 or later
# A copy of GNU AGPL v3 should have been included in this software package in LICENSE.txt.


""" Classifies gene functions according to a curated set of HMM profiles.
    Phylogenetic trees of the gene functions and classification of a gene are
    possible, though must be enabled specifically in the options.
"""

import logging
import os
from typing import Any, Dict, List, Optional

from antismash.common import path
from antismash.common.module_results import ModuleResults
from antismash.common.secmet import Record
from antismash.config import ConfigType
from antismash.config.args import ModuleArgs

from .trees import generate_trees

NAME = "smcogs"
SHORT_DESCRIPTION = "gene function classification via smCOG"


class SMCOGTreeResults(ModuleResults):
    """ Results for the smcogs trees. Tracks the location of the tree images
        generated but does not keep a full copy of the tree.
    """
    schema_version = 2

    def __init__(self, record_id: str, relative_path: str, tree_images: Dict[str, str]) -> None:
        super().__init__(record_id)
        self.tree_images = tree_images  # cds_name -> tree filename
        self.relative_tree_path = relative_path  # path where tree images are saved

    def to_json(self) -> Dict[str, Any]:
        return {"schema_version": self.schema_version,
                "record_id": self.record_id,
                "tree_paths": self.tree_images,
                "image_dir": self.relative_tree_path}

    @staticmethod
    def from_json(json: Dict[str, Any], record: Record) -> Optional["SMCOGTreeResults"]:
        if json.get("schema_version") != SMCOGTreeResults.schema_version:
            logging.debug("Schema version mismatch, discarding SMCOGs results")
            return None
        if record.id != json.get("record_id"):
            logging.debug("Record ID mismatch, discarding SMCOGs results")
            return None
        return SMCOGTreeResults(json["record_id"], json["image_dir"], json["tree_paths"])

    def add_to_record(self, record: Record) -> None:
        """ Annotate smCOG tree paths in CDS features """
        logging.debug("annotating genes with SMCOG trees: %d genes", len(self.tree_images))
        for feature in record.get_cds_features_within_regions():
            if feature.get_name() in self.tree_images:
                feature.notes.append(f"smCOG tree PNG image: smcogs/{self.tree_images[feature.get_name()]}")


def check_options(_options: ConfigType) -> List[str]:
    """ Checks options for problems. """
    return []


def get_arguments() -> ModuleArgs:
    """ Construct the arguments.
        Classification is enabled by default, but an extra option for generating
        trees is also required.
    """
    group = ModuleArgs("Basic analysis options", "smcog", basic_help=True)
    group.add_analysis_toggle('--smcog-trees',
                              dest='smcog_trees',
                              action='store_true',
                              default=False,
                              help="Generate phylogenetic trees of sec. "
                                   "met. cluster orthologous groups.")
    return group


def is_enabled(options: ConfigType) -> bool:
    """ Enabled if tree generation is requested """
    return options.smcog_trees


def regenerate_previous_results(results: Dict[str, Any], record: Record, _options: ConfigType
                                ) -> Optional[SMCOGTreeResults]:
    """ Reconstructs the previous results, unless the trees weren't generated
        previously or a previously generated tree output file is missing.
    """
    if not results:
        return None
    parsed = SMCOGTreeResults.from_json(results, record)
    if not parsed:
        return None

    for tree_filename in parsed.tree_images.values():
        if not os.path.exists(os.path.join(parsed.relative_tree_path, tree_filename)):
            logging.debug("Tree image files missing and must be regenerated")
            return None
    return parsed


def check_prereqs(options: ConfigType) -> List[str]:
    "Check if all required applications are around"
    failure_messages = []
    for binary_name in ['fasttree']:
        if binary_name not in options.executables:
            failure_messages.append(f"Failed to locate executable: {binary_name}")

    return failure_messages


def run_on_record(record: Record, results: Optional[SMCOGTreeResults],
                  options: ConfigType) -> SMCOGTreeResults:
    """ Generates phylogeny trees of the classifications made by SMCOGs
    """
    if results and isinstance(results, SMCOGTreeResults):
        return results
    # create the smcogs output directory if required
    relative_output_dir = os.path.relpath(os.path.join(options.output_dir, "smcogs"), os.getcwd())
    smcogs_dir = os.path.abspath(relative_output_dir)
    if not os.path.exists(smcogs_dir):
        os.mkdir(smcogs_dir)

    nrpspks_genes = record.get_nrps_pks_cds_features()
    with path.changed_directory(smcogs_dir):
        trees = generate_trees(smcogs_dir, record.get_cds_features_within_regions(), nrpspks_genes)

    return SMCOGTreeResults(record.id, relative_output_dir, trees)
