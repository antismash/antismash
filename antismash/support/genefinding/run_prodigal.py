# License: GNU Affero General Public License v3 or later
# A copy of GNU AGPL v3 should have been included in this software package in LICENSE.txt.

""" Gene finding using Prodigal
"""

import logging
from os import path

from Bio.SeqFeature import FeatureLocation
from helperlibs.wrappers.io import TemporaryDirectory

from antismash.common.fasta import write_fasta
from antismash.common.secmet import CDSFeature, Record
from antismash.common.secmet.features.cds_feature import MAX_TRANSLATION_LENGTH
from antismash.common.subprocessing import execute
from antismash.config import ConfigType


def run_prodigal(record: Record, options: ConfigType) -> None:
    """ Run progidal to annotate prokaryotic sequences
    """
    if "basedir" in options.get('prodigal', ''):
        basedir = options.prodigal.basedir
    else:
        basedir = ""
    with TemporaryDirectory(change=True):
        name = record.id.lstrip('-')
        if not name:
            name = "unknown"
        fasta_file = f"{name}.fasta"
        result_file = f"{name}.predict"
        write_fasta([name], [str(record.seq)], fasta_file)

        # run prodigal
        prodigal = [path.join(basedir, 'prodigal')]
        prodigal.extend(['-i', fasta_file, '-f', 'sco', '-o', result_file])
        if options.genefinding_tool == "prodigal-m" or len(record.seq) < 20000:
            prodigal.extend(['-p', 'meta'])

        err = execute(prodigal).stderr
        if err.find('Error') > -1:
            logging.error("Failed to run prodigal: %s", err)
            raise RuntimeError(f"prodigal error: {err}")
        found = 0
        with open(result_file, "r") as handle:
            lines = handle.readlines()
        for line in lines:
            # skip first line
            if not line.startswith('>'):
                continue
            name, start_chunk, end_chunk, prodigal_strand = line[1:].rstrip().split("_")

            try:
                start = int(start_chunk)
                end = int(end_chunk)
                if prodigal_strand == "+":
                    strand = 1
                else:
                    strand = -1
            except ValueError:
                logging.error('Malformatted prodigal output line %r', line.rstrip())
                continue

            if start > end:
                strand = -1
                start, end = end, start

            length = end - start
            if length > (MAX_TRANSLATION_LENGTH * 3) * .9:
                logging.warning("Ignoring potential gene too long for dependencies: %dkb", length // 1000)
                continue

            loc = FeatureLocation(start-1, end, strand=strand)
            translation = record.get_aa_translation_from_location(loc)
            feature = CDSFeature(loc, locus_tag=f"ctg{record.record_index}_{name}",
                                 translation=translation, translation_table=record.transl_table)
            record.add_cds_feature(feature)
            found += 1
    logging.debug("prodigal found %d CDS features", found)
