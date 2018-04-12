# License: GNU Affero General Public License v3 or later
# A copy of GNU AGPL v3 should have been included in this software package in LICENSE.txt.

"""SVG output format module

"""
import logging
import os
from typing import Dict, List

from antismash.common.module_results import ModuleResults
from antismash.config import ConfigType
from antismash.modules import clusterblast

NAME = "svg"
SHORT_DESCRIPTION = "SVG output"


def write(options: ConfigType, results: List[Dict[str, ModuleResults]]) -> None:
    """ Writes all SVG files into the output directory """
    svgdir = os.path.join(options.output_dir, "svg")
    logging.debug("Writing seq_records SVGs to %r", svgdir)
    if not os.path.exists(svgdir):
        os.mkdir(svgdir)
    for record_result in results:
        result = record_result.get(clusterblast.__name__)
        if not result:
            continue
        assert isinstance(result, clusterblast.ClusterBlastResults), type(result)
        result.write_svg_files(svgdir)
