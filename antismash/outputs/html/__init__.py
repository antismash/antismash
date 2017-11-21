# License: GNU Affero General Public License v3 or later
# A copy of GNU AGPL v3 should have been included in this software package in LICENSE.txt.

"""HTML output format module

"""
import logging
import os
import shutil

from argparse import Namespace

from antismash.common import path, deprecated
from antismash.config.args import ModuleArgs
from antismash.outputs.html.generator import generate_webpage

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
    generate_webpage(records, results, options)


def copy_template_dir(template, output_dir):
    """ Copy files from a template directory to the output directory, removes
        any existing directory first
    """
    target_dir = os.path.join(output_dir, template)
    if os.path.exists(target_dir):
        shutil.rmtree(target_dir)
    shutil.copytree(path.get_full_path(__file__, template), target_dir)
