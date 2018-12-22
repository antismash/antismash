# License: GNU Affero General Public License v3 or later
# A copy of GNU AGPL v3 should have been included in this software package in LICENSE.txt.

"""HTML output format module

"""
import os
import shutil
from typing import Dict, List

from antismash.common import path
from antismash.common.module_results import ModuleResults
from antismash.common.secmet import Record
from antismash.config import ConfigType
from antismash.config.args import ModuleArgs
from antismash.outputs.html.generator import generate_webpage

NAME = "html"
SHORT_DESCRIPTION = "HTML output"


def get_arguments() -> ModuleArgs:
    """ Builds the arguments for the HMTL output module """
    args = ModuleArgs("Output options", "html", enabled_by_default=True)
    args.add_option("--html-title",
                    dest="html_title",
                    type=str,
                    default="",
                    help=("Custom title for the HTML output page "
                          "(default is input filename)."))
    args.add_option("--html-description",
                    dest="html_description",
                    type=str,
                    default="",
                    help="Custom description to add to the output.")
    return args


def check_options(_options: ConfigType) -> List[str]:
    """ Check options, but none to check """
    return []


def is_enable(_options: ConfigType) -> bool:
    """ Is the HMTL module enabled (currently always enabled) """
    return True  # TODO: add an arg to disable


def write(records: List[Record], results: List[Dict[str, ModuleResults]],
          options: ConfigType) -> None:
    """ Writes all results to a webpage, where applicable. Writes to options.output_dir

        Arguments:
            records: the list of Records for which results exist
            results: a list of dictionaries containing all module results for records
            options: antismash config object

        Returns:
            None
    """
    output_dir = options.output_dir

    copy_template_dir('css', output_dir)
    copy_template_dir('js', output_dir)
    copy_template_dir('images', output_dir)

    # Generate structure images for records obtained from BioSQL
    generate_webpage(records, results, options)


def copy_template_dir(template: str, output_dir: str) -> None:
    """ Copy files from a template directory to the output directory, removes
        any existing directory first
    """
    target_dir = os.path.join(output_dir, template)
    if os.path.exists(target_dir):
        shutil.rmtree(target_dir)
    shutil.copytree(path.get_full_path(__file__, template), target_dir)
