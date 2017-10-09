# License: GNU Affero General Public License v3 or later
# A copy of GNU AGPL v3 should have been included in this software package in LICENSE.txt.

import logging

from antismash.common.all_orfs import find_all_orfs

from .run_prodigal import run_prodigal
from .run_glimmerhmm import run_glimmerhmm


def run_on_record(record, options):
    "Find genes in a seq_record"
    if options.genefinding_tool == 'none':
        raise ValueError("Called find_genes, but genefinding disabled")
    if options.genefinding_tool == "all-orfs":
        find_all_orfs(record)
    elif options.taxon == 'fungi':
        logging.debug("Running glimmerhmm genefinding")
        run_glimmerhmm(record, options)
    elif options.genefinding_tool in ["prodigal", "prodigal-m"]:
        logging.debug("Running prodigal based genefinding")
        run_prodigal(record, options)
