# License: GNU Affero General Public License v3 or later
# A copy of GNU AGPL v3 should have been included in this software package in LICENSE.txt.

"""SVG output format module

"""
import logging
import os
import shutil

import antismash.common.deprecated as utils

NAME = "svg"
SHORT_DESCRIPTION = "SVG output"
ENABLED = True

def write(seq_records, options, results):
    svgdir = os.path.join(options.output_dir, "svg")
    logging.debug("Writing seq_records SVGs to %r", svgdir)
    if not os.path.exists(svgdir):
        os.mkdir(svgdir)
    for seq_record, record_result in zip(seq_records, results):
        result = record_result["modules"].get("antismash.modules.clusterblast")
        if not result:
            continue
        result.write_svg_files(svgdir)

