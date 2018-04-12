# License: GNU Affero General Public License v3 or later
# A copy of GNU AGPL v3 should have been included in this software package in LICENSE.txt.

""" MEME and FIMO runners which don't belong in common.subprocessing """

import os

from antismash.common.subprocessing import parallel_execute
from antismash.common.secmet import Record
from antismash.config import ConfigType


def run_meme(meme_dir: str, _options: ConfigType, verbose: bool) -> int:
    """Set paths, check existing files and run MEME in parallel on each promoter set"""
    args = []
    for plus_minus in os.listdir(meme_dir):
        input_file = os.path.join(meme_dir, plus_minus, "promoters.fasta")
        output_file = os.path.join(meme_dir, plus_minus, "meme.xml")

        # input file present and size not zero --> should be fine to use
        if os.path.isfile(input_file) and os.path.getsize(input_file) > 0:
            # output file already present and size not zero --> do not run MEME again on this one
            if os.path.isfile(output_file) and os.path.getsize(output_file) > 0:
                pass
            else:
                args.append([
                    "meme", input_file,
                    "-oc", os.path.dirname(input_file),
                    "-dna",
                    "-nostatus",
                    "-mod", "anr",
                    "-nmotifs", "1",
                    "-minw", "6",
                    "-maxw", "12",
                    "-revcomp",
                    "-evt", "1.0e+005",
                ])

    errors = parallel_execute(args, verbose=verbose)
    return sum(errors)


def run_fimo(meme_dir: str, fimo_dir: str, record: Record, options: ConfigType, verbose: bool) -> int:
    """Set paths, check existing files and run FIMO in parallel on each predicted motif"""
    if not os.path.exists(fimo_dir):
        os.makedirs(fimo_dir)

    args = []
    for plus_minus in os.listdir(meme_dir):
        motif_file = os.path.join(meme_dir, plus_minus, "meme.html")
        sites_file = os.path.join(meme_dir, plus_minus, "binding_sites.fasta")
        output_file = os.path.join(fimo_dir, plus_minus, "fimo.txt")
        output_dir = os.path.join(fimo_dir, plus_minus)

        # input file and binding sites file present and size not zero --> should be fine to use
        if (os.path.isfile(motif_file)
                and os.path.getsize(motif_file) > 0
                and os.path.isfile(sites_file)
                and os.path.getsize(sites_file) > 0):
            # output file already present and size not zero --> do not run FIMO again on this one
            if os.path.isfile(output_file) and os.path.getsize(output_file) > 0:
                pass
            else:
                args.append([
                    "fimo",
                    "-verbosity", "1",
                    "-motif", "1",
                    "-thresh", "0.00006",
                    "-oc", output_dir,
                    motif_file,
                    os.path.join(options.output_dir, record.name + "_promoter_sequences.fasta"),
                ])

    errors = parallel_execute(args, verbose=verbose)
    return sum(errors)
