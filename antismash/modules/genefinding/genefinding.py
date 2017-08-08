# License: GNU Affero General Public License v3 or later
# A copy of GNU AGPL v3 should have been included in this software package in LICENSE.txt.

from .all_orfs import find_all_orfs
from .run_prodigal import run_prodigal
from .run_glimmer import run_glimmer
from .run_glimmerhmm import run_glimmerhmm

def run_on_record(record, options):
    "Find genes in a seq_record"
    if options.genefinding_tool == 'none':
        raise ValueError("Called find_genes, but genefinding disabled")
    if options.genefinding_all_orfs:
        find_all_orfs(record)
    if options.taxon == 'fungi':
        run_glimmerhmm(record, options)
        return
    if options.genefinding_tool in ["prodigal", "prodigal-m"]:
        run_prodigal(record, options)
    else:
        run_glimmer(record, options)
