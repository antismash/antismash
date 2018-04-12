# License: GNU Affero General Public License v3 or later
# A copy of GNU AGPL v3 should have been included in this software package in LICENSE.txt.

"""Gene finding using GlimmerHMM

   mostly for fungi/eukaryotes
"""

from io import StringIO
import logging

from helperlibs.bio import seqio
from helperlibs.wrappers.io import TemporaryDirectory

from antismash.common import path
from antismash.common.gff_parser import get_features_from_file
from antismash.common.secmet import Record
from antismash.common.subprocessing import execute


def write_search_fasta(record: Record) -> str:
    """ Constructs a FASTA representation of a record and writes it to a
        file in the current directory.

        Returns:
            the name of the file created
    """
    name = record.id.lstrip('-')
    if not name:
        name = "unknown"
    filename = "{}.fasta".format(name)
    with open(filename, 'w') as handle:
        seqio.write([record.to_biopython()], handle, 'fasta')
    return filename


def run_external(fasta_filename: str) -> str:
    """ Runs glimmerhmm on the provided fasta file and returns the stdout output
        from glimmerhmm.
    """
    glimmerhmm = ['glimmerhmm', fasta_filename,
                  path.get_full_path(__file__, "data/train_crypto"), "-g"]
    run_result = execute(glimmerhmm)
    if run_result.stderr.find('ERROR') > -1:
        logging.error("Failed to run GlimmerHMM: %r", run_result.stderr)
        raise RuntimeError("Failed to run GlimmerHMM: %s" % run_result.stderr)
    if "CDS" not in run_result.stdout:
        logging.error("GlimmerHMM gene prediction failed: no genes found.")
        raise RuntimeError("GlimmerHMM found no genes")
    return run_result.stdout


def run_glimmerhmm(record: Record) -> None:
    """ Run glimmerhmm on the record, parse the results and add all detected
        genes to the record
    """
    with TemporaryDirectory(change=True):
        # Write FASTA file and run GlimmerHMM
        fasta_file = write_search_fasta(record)
        results_text = run_external(fasta_file)

    handle = StringIO(results_text)
    features = get_features_from_file(record, handle)
    for feature in features:
        record.add_biopython_feature(feature)
