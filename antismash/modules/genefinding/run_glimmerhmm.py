# vim: set fileencoding=utf-8 :
#
# Copyright (C) 2010-2012 Marnix H. Medema
# University of Groningen
# Department of Microbial Physiology / Groningen Bioinformatics Centre
#
# Copyright (C) 2011,2012 Kai Blin
# University of Tuebingen
# Interfaculty Institute of Microbiology and Infection Medicine
# Div. of Microbiology/Biotechnology
#
# License: GNU Affero General Public License v3 or later
# A copy of GNU AGPL v3 should have been included in this software package in LICENSE.txt.

"""Gene finding using GlimmerHMM


mostly for fungi/eukaryotes
"""

import logging
from helperlibs.wrappers.io import TemporaryDirectory
from helperlibs.bio import seqio
from Bio.SeqFeature import SeqFeature, FeatureLocation, CompoundLocation

from antismash.common import deprecated as utils
from antismash.common import path
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
    with TemporaryDirectory(change=True):

        utils.fix_record_name_id(seq_record, options)
        #Write FASTA file and run GlimmerHMM
        fasta_file = write_search_fasta(seq_record)
        resultstext = run_external(fasta_file)

        #Parse GlimmerHMM predictions
        resultstext = resultstext.replace("\r", " ")
        lines = resultstext.split("\n")[2:-1]
        orfnames = []
        positions = []
        strands = []
        x = 0
        orfnr = 0
        starts = []
        ends = []
        for x, line in enumerate(lines):
            columns = line.split("\t")
            if len(columns) < 7:
                raise ValueError("GlimmerHMM result line too short: %s" % line)

            if x == 0:
                if columns[6] == "+":
                    bpy_strand = 1
                else:
                    bpy_strand = -1
                if "mRNA" not in line:
                    starts.append(int(columns[3]))
                    ends.append(int(columns[4]))
            elif x == (len(lines) - 1) or "mRNA" in lines[x + 1]:
                if columns[6] == "+":
                    bpy_strand = 1
                else:
                    bpy_strand = -1
                strands.append(bpy_strand)
                starts.append(int(columns[3]))
                ends.append(int(columns[4]))
                orfnames.append("orf" + (5 - orfnr) * "0" + str(orfnr))
                orfnr += 1
                if len(starts) == 1:
                    if starts[0] == 0:
                        starts[0] = 1
                    if ends[0] == 0:
                        ends[0] = 1
                    positions.append([[starts[0] - 1, ends[0]]])
                else:
                    pos = []
                    if bpy_strand == -1:
                        starts.reverse()
                        ends.reverse()
                    for i in starts:
                        if i == 0:
                            i = 1
                        if ends[starts.index(i)] == 0:
                            ends[starts.index(i)] = 1
                        pos.append([i - 1, ends[starts.index(i)]])
                    positions.append(pos)
                starts = []
                ends = []
            elif "mRNA" not in line:
                starts.append(int(columns[3]))
                ends.append(int(columns[4]))
        if not orfnames:
            logging.error("GlimmerHMM gene prediction failed. Please check the " \
                "format of your input FASTA file.")
        #Create seq_record features for identified genes
        idx = 0
        for orfname in orfnames:
            bpy_strand = strands[idx]
            genepositions = positions[idx]
            #For genes with only one CDS
            if len(genepositions) == 1:
                gstart, gend = genepositions[0]
                loc = FeatureLocation(gstart, gend, strand=bpy_strand)
                feature = SeqFeature(location=loc, id=orfname, type="CDS",
                            qualifiers={'locus_tag': ['ctg%s_%s' % (options.record_idx, orfname)]})
                seq_record.features.append(feature)
            #For genes with multiple exons
            else:
                gstart, gend = min(genepositions[0]), max(genepositions[-1])
                sublocations = []
                for exonstart, exonend in genepositions:
                    exonloc = FeatureLocation(exonstart, exonend, strand=bpy_strand)
                    sublocations.append(exonloc)
                loc = CompoundLocation(sublocations)
                feature = SeqFeature(location=loc, id=orfname, type="CDS",
                            qualifiers={'locus_tag': ['ctg%s_%s' % (options.record_idx, orfname)]})
                seq_record.features.append(feature)
            idx += 1
