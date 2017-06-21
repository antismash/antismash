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

"""ORF detector for all ORFs > 60 nt

"""

import logging
from Bio.SeqFeature import SeqFeature, FeatureLocation

class Orf:
    """A putative open reading frame"""
    def __init__(self, start, stop, direction):
        self.start = start
        self.stop = stop
        self.direction = direction

    def __str__(self):
        dir_string = "+"
        if self.direction == -1:
            dir_string = "-"
        return "%9d%9d  %s%d" % (self.start+1, self.stop, dir_string, self.get_frame())

# TODO: python3ize
#    def __cmp__(self, other):
#        return cmp(self.start, other.start)

    def get_frame(self):
        """Calculate the ORF's frame offset"""
        return (self.start % 3) + 1

    def len(self):
        """Get the length of the Orf"""
        return self.stop+2 - self.start

def scan_orfs(seq, direction, offset=0):
    """Scan for open reading frames on a given sequence"""
    seq = seq.upper()
    START_CODONS = ('ATG', 'GTG', 'TTG')
    STOP_CODONS  = ('TAA', 'TAG', 'TGA')
    matches = []
    # Remember the last stop codon found per frame, so we can take some
    # shortcuts later
    last_stop = [0, 0, 0]
    last_orf = [None, None, None]
    # cache the sequence length
    seq_len = len(seq)
    for i in range(0, seq_len - 2 ):
        if seq[i:i+3] in START_CODONS or i == 0:
            # If the last stop codon found is in the frame of this start codon
            # and the start codon is upstream of the stop codon, we have
            # already found a start codon further upstream.
            if i < last_stop[i%3]:
                continue
            # Look for the next stop codon in this frame
            do_continue = False
            for j in range(i, seq_len - 2, 3):
                if seq[j:j+3] in STOP_CODONS:
                    # Skip Orfs that are shorter than 20 AA / 60 bases
                    if j - i > 60:
                        if direction == 1:
                            new_orf = Orf(i + offset, j + 2 + offset, direction)
                        else:
                            # i and j are the position on the reverse strand, convert this back to
                            # the forward strand positions
                            new_orf = Orf(seq_len - (j + 2) + offset, seq_len - i + offset, direction)
                        matches.append(new_orf)
                        # This was a good hit, update the last_stop cache.
                        last_stop[j % 3] = j
                        last_orf[i % 3] = new_orf
                    do_continue = True
                    break
            if do_continue:
                continue
            #Save orfs ending at the end of the sequence without stop codon
            j = seq_len - 1
            if direction == 1:
                new_orf = Orf(i + offset, j + 2 + offset, direction)
            else:
                # i and j are the position on the reverse strand, convert this back to
                # the forward strand positions
                new_orf = Orf(seq_len - (j + 2) + offset, seq_len - i + offset, direction)
            matches.append(new_orf)
            last_stop[i % 3] = j+2
            last_orf[i % 3] = new_orf
    return matches

def get_reverse_complement(seq):
    """generate the reverse strand of a given sequence"""
    seq = seq.upper()
    reverse_seq = ""
    for i in range(len(seq) - 1, -1, -1):
        if seq[i] == 'G':
            reverse_seq += 'C'
        elif seq[i] == 'C':
            reverse_seq += 'G'
        elif seq[i] == 'A':
            reverse_seq += 'T'
        elif seq[i] == 'T':
            reverse_seq += 'A'
        else:
            reverse_seq += 'X'
    return reverse_seq

def sort_orfs(orfs):
    if len(orfs) == 0:
        return orfs
    startpositions = [min([orf.start, orf.stop]) for orf in orfs]
    positions_and_orfs = sorted(zip(startpositions, orfs))
    startpositions, orfs = zip(*positions_and_orfs)
    return orfs

def find_all_orfs(seq_record, options):
    #Get sequence
    fasta_seq = str(seq_record.seq)
    #Find orfs
    forward_matches = scan_orfs(fasta_seq, 1)
    reverse_matches = scan_orfs(get_reverse_complement(fasta_seq), -1)
    all_orfs = forward_matches + reverse_matches
    #Create seq_record features for identified genes
    if len(all_orfs) == 0:
        logging.error("No ORFs found. Please check the " \
            "format of your input FASTA file.")
    orfnr = 0
    for orf in sort_orfs(all_orfs):
        if orf.start < 0:
            orf.start = 0
        if orf.stop < 0:
            orf.stop = 0
        seqlength = len(str(seq_record.seq))
        while orf.start > seqlength:
            orf.start = orf.start - 3
        while orf.stop > seqlength:
            orf.stop = orf.stop - 3
        loc = FeatureLocation(orf.start, orf.stop, strand=orf.direction)
        feature = SeqFeature(location=loc, id=str(orf), type="CDS",
                    qualifiers={'locus_tag': ['ctg%s_allorf%s%s' % (options.record_idx, "0" * (6 - len(str(orfnr))), str(orfnr))]})
        feature.qualifiers['note'] = ["auto-all-orf"]
        seq_record.features.append(feature)
        orfnr += 1
