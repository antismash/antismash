# License: GNU Affero General Public License v3 or later
# A copy of GNU AGPL v3 should have been included in this software package in LICENSE.txt.

"""Gene finding using GlimmerHMM

   mostly for fungi/eukaryotes
"""

import logging
from helperlibs.wrappers.io import TemporaryDirectory
from helperlibs.bio import seqio
from io import StringIO
from Bio.SeqFeature import SeqFeature, FeatureLocation, CompoundLocation

from antismash.common import deprecated as utils
from antismash.common import path
from antismash.common.gff_parser import get_features_from_file
from antismash.common.subprocessing import execute


def write_search_fasta(seq_record):
    name = seq_record.id.lstrip('-')
    if not name:
        name = "unknown"
    filename = "{}.fasta".format(name)
    with open(filename, 'w') as handle:
        seqio.write([seq_record], handle, 'fasta')
    return filename

def run_external(fasta_filename):
    glimmerhmm = ['glimmerhmm', fasta_filename,
                  path.get_full_path(__file__, "data/train_crypto"), "-g"]
    run_result = execute(glimmerhmm)
    if run_result.stderr.find('ERROR') > -1:
        logging.error("Failed to run GlimmerHMM: %r", run_result.stderr)
        raise RuntimeError("Failed to run GlimmerHMM: %s", run_result.stderr)
    if "CDS" not in run_result.stdout:
        logging.error("GlimmerHMM gene prediction failed: no genes found.")
        raise RuntimeError("GlimmerHMM found no genes")
    return run_result.stdout

def run_glimmerhmm(seq_record, options):
    utils.fix_record_name_id(seq_record, options)
    with TemporaryDirectory(change=True):
        # Write FASTA file and run GlimmerHMM
        fasta_file = write_search_fasta(seq_record)
        results_text = run_external(fasta_file)

    handle = StringIO(results_text)
    features = get_features_from_file(seq_record, handle)
    seq_record.features.extend(features)
