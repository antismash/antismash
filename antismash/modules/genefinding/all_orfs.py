# License: GNU Affero General Public License v3 or later
# A copy of GNU AGPL v3 should have been included in this software package in LICENSE.txt.

import logging
from Bio.SeqFeature import SeqFeature, FeatureLocation, \
                           BeforePosition, AfterPosition

def scan_orfs(seq, direction, offset=0):
    """ Scan for open reading frames on a given sequence
        skips all ORFs with a size less than 60 bases
    """
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
            matches.append(new_orf)
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

        # if we found a matching stop, carry on looking for starts
        if last_stop[i % 3] > i:
            continue

        # Save orfs ending at the end of the sequence without stop codon
        if direction == 1:
            new_orf = FeatureLocation(i + offset, AfterPosition(seq_len + offset), direction)
        else:
            # reversed, so convert back to the forward positions
            new_orf = FeatureLocation(BeforePosition(offset), offset + seq_len - i, direction)
        matches.append(new_orf)
        # since there are no stop codons, just stop here
        break
    return matches

def sort_orfs(orfs):
    startpositions = [min([orf.start, orf.end]) for orf in orfs]
    positions_and_orfs = sorted(zip(startpositions, orfs), key=lambda x: x[0])
    startpositions, orfs = zip(*positions_and_orfs)
    return orfs

def find_all_orfs(seq_record):
    logging.debug("Finding all ORFs")
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
        return len(all_orfs)
    orfs = sort_orfs(all_orfs)
    for orfnr, orf in enumerate(orfs):
        locus_tag = 'ctg%d_allorf%06d' % (seq_record.record_index, orfnr)
        feature = SeqFeature(location=orf, id=locus_tag, type="CDS",
                    qualifiers={'locus_tag': [locus_tag]})
        feature.qualifiers['note'] = ["auto-all-orf"]
        seq_record.features.append(feature)
    logging.info("Found %d ORFs", len(orfs))
    return len(orfs)
