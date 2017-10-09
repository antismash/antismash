# License: GNU Affero General Public License v3 or later
# A copy of GNU AGPL v3 should have been included in this software package in LICENSE.txt.

"""HTML output format module

"""
import logging
from os import path
import shutil

from antismash.outputs.html.generator import generate_webpage
from antismash.config.args import ModuleArgs
from antismash.outputs.html.structure_drawer import generate_chemical_structure_preds

NAME = "html"
SHORT_DESCRIPTION = "HTML output"


def get_arguments():
    return ModuleArgs("Output options", "html", enabled_by_default=True)


def check_options(options):
    return []


def is_enabled(options):
    return True  # TODO: add an arg to disable


def write(records, results, options):
    output_dir = options.output_dir

    copy_template_dir('css', output_dir)
    copy_template_dir('js', output_dir)
    copy_template_dir('images', output_dir)

    # Generate structure images for records obtained from BioSQL
    generate_structure_images(records, options)
    generate_webpage(records, results, options)


def copy_template_dir(template, output_dir):
    "Copy files from a template directory to the output directory"
    basedir = path.dirname(__file__)

    target_dir = path.join(output_dir, template)
    if path.exists(target_dir):
        shutil.rmtree(target_dir)
    shutil.copytree(path.join(basedir, template), target_dir)


def generate_structure_images(records, options):
    "Generate the structure images based on Monomers prediction in cluster feature"

    logging.critical("pksnrps results would be added here, but shouldn't be")
    return

#    for record in records:
#        # Ugly temporary solution:
#        # At first we have to regenerate the relevant information for the pksnrpsvars dictionary from the record file
#        pksnrpsvars = utils.Storage()
#        pksnrpsvars.compound_pred_dict = {}
#        pksnrpsvars.failedstructures = []

#        for cluster in record.get_clusters():
#            cluster_number = cluster.get_cluster_number()
#            pksnrpsvars.compound_pred_dict[cluster_number] = utils.get_structure_pred(cluster)
#        if pksnrpsvars.compound_pred_dict:
#            generate_chemical_structure_preds(pksnrpsvars, record, options)
