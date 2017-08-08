# License: GNU Affero General Public License v3 or later
# A copy of GNU AGPL v3 should have been included in this software package in LICENSE.txt.

"""Gene finding using Glimmer

"""

import logging
from os import path
from Bio.SeqFeature import SeqFeature, FeatureLocation
from helperlibs.wrappers.io import TemporaryDirectory
from helperlibs.bio import seqio

from antismash.common import deprecated as utils
from antismash.common.subprocessing import execute

def run_glimmer(seq_record, options):
    "Run glimmer3 to annotate prokaryotic sequences"
    basedir = options.get('glimmer', {}).get('basedir', '')
    with TemporaryDirectory(change=True):
        name = seq_record.id.lstrip('-')
        if not name:
            name = "unknown"
        fasta_file = '%s.fasta' % name
        longorfs_file = '%s.longorfs' % name
        icm_file = '%s.icm' % name
        result_file = '%s.predict' % name

        # run long-orfs
        with open(fasta_file, 'w') as handle:
            seqio.write([seq_record], handle, 'fasta')
        long_orfs = [path.join(basedir, 'long-orfs')]
        long_orfs.extend(['-l', '-n', '-t', '1.15', '--trans_table', '11',
                          fasta_file, longorfs_file])
        run_result = execute(long_orfs)
        if run_result.stderr.find('ERROR') > -1:
            logging.error("Locating long orfs failed: %r", run_result.stderr)
            return

        # run extract
        extract = [path.join(basedir, 'extract'), '-t', fasta_file,
                   longorfs_file]
        run_result = execute(extract)
        if not run_result.stdout:
            logging.error("Failed to extract genes from model, aborting: %r", run_result.stderr)
            return

        build_icm = [path.join(basedir, 'build-icm'), '-r', icm_file]
        run_result = execute(build_icm, stdin=run_result.stdout)
        if not run_result.stderr:
            logging.error("Failed to build gene model: %r", run_result.stderr)
            return

        # run glimmer3
        glimmer = [path.join(basedir, 'glimmer3')]
        glimmer.extend(['-l', '-o', '50', '-g', '90', '-q', '3000', '-t', '30',
                        '--trans_table', '11', fasta_file, icm_file, name])

        run_result = execute(glimmer)
        if run_result.stderr.find('ERROR') > -1:
            logging.error("Failed to run glimmer3: %r", run_result.stderr)
            return
        for line in open(result_file, 'r'):
            # skip first line
            if line.startswith('>'):
                continue

            name, start, end, strand, score = line.split()

            try:
                start = int(start)
                end = int(end)
                strand = int(strand)
            except ValueError:
                logging.error('Malformatted glimmer output line %r', line.rstrip())

            if start > end:
                bpy_strand = -1
                tmp = start
                start = end
                end = tmp
            else:
                bpy_strand = 1

            loc = FeatureLocation(start-1, end, strand=bpy_strand)
            feature = SeqFeature(location=loc, id=name, type="CDS",
                        qualifiers={'locus_tag': ['ctg%s_%s' % (seq_record.record_index, name)],
                                    'note': ['Glimmer score: %s' %score]})
            seq_record.features.append(feature)
