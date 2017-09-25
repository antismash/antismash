# License: GNU Affero General Public License v3 or later
# A copy of GNU AGPL v3 should have been included in this software package in LICENSE.txt.

"""HTML output format module

"""
import logging
from os import path
import shutil

from antismash.outputs.html.generator import generate_webpage
import antismash.common.deprecated as utils
from antismash.config.args import ModuleArgs
from antismash.outputs.html.structure_drawer import generate_chemical_structure_preds

NAME = "html"
SHORT_DESCRIPTION = "HTML output"

def get_arguments():
    return ModuleArgs("Output options", "html", enabled_by_default=True)

def check_options(options):
    return []

def is_enabled(options):
    return True #TODO: add an arg to disable

def write(seq_records, results, options):
    output_dir = options.output_dir

    copy_template_dir('css', output_dir)
    copy_template_dir('js', output_dir)
    copy_template_dir('images', output_dir)

    # Generate structure images for records obtained from BioSQL
    generate_structure_images(seq_records, options)
    generate_webpage(seq_records, results, options)

def copy_template_dir(template, output_dir):
    "Copy files from a template directory to the output directory"
    basedir = path.dirname(__file__)

    target_dir = path.join(output_dir, template)
    if path.exists(target_dir):
        shutil.rmtree(target_dir)
    shutil.copytree(path.join(basedir, template), target_dir)

def generate_structure_images(seq_records, options):
    "Generate the structure images based on Monomers prediction in cluster feature"

    logging.critical("pksnrps results would be added here, but shouldn't be")
    return

#    for seq_record in seq_records:
#        # Ugly temporary solution:
#        # At first we have to regenerate the relevant information for the pksnrpsvars dictionary from the seq_record file
#        pksnrpsvars = utils.Storage()
#        pksnrpsvars.compound_pred_dict = {}
#        pksnrpsvars.failedstructures = []

#        geneclusters = utils.get_cluster_features(seq_record)

#        for genecluster in geneclusters:
#            geneclusternr = utils.get_cluster_number(genecluster)
#            pksnrpsvars.compound_pred_dict[geneclusternr] = utils.get_structure_pred(genecluster)
#        if len(pksnrpsvars.compound_pred_dict) > 0:
#            generate_chemical_structure_preds(pksnrpsvars, seq_record, options)
