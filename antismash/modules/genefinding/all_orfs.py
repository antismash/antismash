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
import math
from Bio.SeqFeature import SeqFeature, FeatureLocation, \
                           BeforePosition, AfterPosition

def scan_orfs(seq, direction, offset=0):
    """Scan for open reading frames on a given sequence"""
    seq = seq.upper()
    start_codons = ('ATG', 'GTG', 'TTG')
    stop_codons = ('TAA', 'TAG', 'TGA')
    matches = []
    # Remember the last stop codon found per frame, so we can take some
    # shortcuts later
    last_stop = [0, 0, 0]
    # cache the sequence length
    seq_len = len(seq)
    for i in range(0, seq_len - 2):
        # If the last stop codon found is in the frame of this start codon
        # and the start codon is upstream of the stop codon, we have
        # already found a start codon further upstream.
        if i < last_stop[i % 3]:
            continue
        if seq[i:i+3] in stop_codons and last_stop[i % 3] == 0:
            # special case for unstarted stops
            last_stop[i % 3] = i
            if direction == 1:
                new_orf = FeatureLocation(BeforePosition(offset), offset + i + 2, direction)
            else:
                new_orf = FeatureLocation(seq_len + offset - (i + 2),
                              AfterPosition(seq_len + offset), direction)
        if seq[i:i+3] not in start_codons:
            continue
        # Look for the next stop codon in this frame
        for j in range(i, seq_len - 2, 3):
            if seq[j:j+3] in stop_codons:
                last_stop[j % 3] = j
                # Skip Orfs that are shorter than 20 AA / 60 bases
                if j - i <= 60:
                    break # since no ORFs will be bigger before the stop
                start = i
                end = j + 2
                if direction == 1:
                    new_orf = FeatureLocation(offset + start,
                                              offset + end, direction)
                else:
                    # reversed, so convert back to the forward positions
                    new_orf = FeatureLocation(seq_len + offset - end,
                                              seq_len + offset - start, direction)
                matches.append(new_orf)
                # This was a good hit, update the last_stop cache.
                break

        # walked off the end without finding a stop
        if j < seq_len - 2:
            continue

        # Save orfs ending at the end of the sequence without stop codon
        end = seq_len + 1
        if direction == 1:
            new_orf = FeatureLocation(i + offset, AfterPosition(end + offset), direction)
        else:
            # reversed, so convert back to the forward positions
            new_orf = FeatureLocation(BeforePosition(offset), offset + seq_len - i, direction)
        matches.append(new_orf)
        # since there are no stop codons, just stop here
        break
    return matches

def sort_orfs(orfs):
    if not orfs:
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
    reverse_matches = scan_orfs(str(seq_record.seq.complement()), -1)
    all_orfs = forward_matches + reverse_matches
    #Create seq_record features for identified genes
    if not all_orfs:
        logging.error("No ORFs found. Please check the " \
            "format of your input FASTA file.")
    orfnr = 0
    for orf in sort_orfs(all_orfs):
        # ensure the location isn't negative to start with
        if orf.start < 0:
            orf.start = 0
        if orf.stop < 0:
            orf.stop = 0
        seqlength = len(str(seq_record.seq))
        # make sure we don't create a negative location if we adjust them
        assert seqlength >= 3 or orf.start < 3 and orf.stop < 3
        if orf.start > seqlength:
            orf.start -= math.ceil((orf.start - seqlength) / 3) * 3
        if orf.stop > seqlength:
            orf.stop -= math.ceil((orf.stop - seqlength) / 3) * 3
        loc = FeatureLocation(orf.start, orf.stop, strand=orf.direction)
        locus_tag = 'ctg%s_allorf%06d' % (options.record_idx, orfnr)
        feature = SeqFeature(location=loc, id=str(orf), type="CDS",
                    qualifiers={'locus_tag': [locus_tag]})
        feature.qualifiers['note'] = ["auto-all-orf"]
        seq_record.features.append(feature)
        orfnr += 1
