# License: GNU Affero General Public License v3 or later
# A copy of GNU AGPL v3 should have been included in this software package in LICENSE.txt.

from collections import OrderedDict
import logging
from typing import Dict, List, Set, Tuple

import Bio.Data.IUPACData
from Bio.SeqUtils.ProtParam import ProteinAnalysis


def read_fasta(filename: str) -> Dict[str, str]:
    """ reads a fasta file into a dict: id -> sequence, returns the dict """
    ids = []
    sequence_info = []
    with open(filename, "r") as fasta:
        current_seq = []
        for line in fasta:
            line = line.strip()
            if not line:
                continue
            if line[0] == '>':
                ids.append(line[1:].replace(" ", "_"))
                if current_seq:
                    sequence_info.append("".join(current_seq))
                    current_seq = []
            else:
                if not ids:
                    raise ValueError("Sequence before identifier in fasta file")
                if not line.replace("-", "z").isalpha():
                    raise ValueError("Sequence contains non-alphabetic characters")
                current_seq.append(line)
    if current_seq:
        sequence_info.append("".join(current_seq))
    if len(ids) != len(sequence_info):
        raise ValueError("Fasta files contains different counts of sequences and ids")
    if not ids:
        logging.debug("Fasta file %s contains no sequences", filename)
        raise ValueError("Fasta file contains no sequences")
    return OrderedDict(zip(ids, sequence_info))


def generate_unique_id(prefix: str, existing_ids: Set, start=0, max_length=-1) -> Tuple[str, int]:
    """ Generate a identifier of the form prefix_num, e.g. seq_15.

        Does *not* add the generated prefix to current identifiers

        Args:
            prefix: The text portion of the name.
            existing_ids: The current identifiers to avoid collision with.
            start: An integer to start counting at (default: 0)
            max_length: The maximum length allowed for the identifier,
                        values less than 1 are considerd to be no limit.

        Returns:
            A tuple of the identifier generated and the value of the counter
                at the time the identifier was generated, e.g. ("seq_15", 15)

    """
    counter = int(start)
    existing_ids = set(existing_ids)
    max_length = int(max_length)

    format_string = "{}_%d".format(prefix)
    name = format_string % counter
    while name in existing_ids:
        counter += 1
        name = format_string % counter
    if max_length > 0 and len(name) > max_length:
        raise RuntimeError("Could not generate unique id for %s after %d iterations" % (prefix, counter - start))
    return name, counter


class RobustProteinAnalysis(ProteinAnalysis):
    """ A simple subclass of ProteinAnalysis that can deal with
        a protein sequence containing invalid characters.

        If ignoring invalid characters, the molecular weight is increased by
        the average weight of an amino-acid (i.e. 110) for each invalid case.
    """
    PROTEIN_LETTERS = set(Bio.Data.IUPACData.protein_letters)

    def __init__(self, prot_sequence, monoisotopic=False, ignore_invalid=True) -> None:
        if not isinstance(ignore_invalid, bool):
            raise TypeError("ignore_invalid must be a boolean")
        self._ignore_invalid = ignore_invalid

        prot_sequence = prot_sequence.upper()

        self.original_sequence = prot_sequence
        # remove all invalids
        prot_sequence = []
        for i in self.original_sequence:
            if i in RobustProteinAnalysis.PROTEIN_LETTERS:
                prot_sequence.append(i)
        prot_sequence = "".join(prot_sequence)
        super(RobustProteinAnalysis, self).__init__(prot_sequence, monoisotopic)

    def molecular_weight(self) -> float:
        weight = super(RobustProteinAnalysis, self).molecular_weight()
        if not self._ignore_invalid:
            aa_difference = len(self.original_sequence) - len(self.sequence)
            weight += 110 * aa_difference
        return weight

def extract_by_reference_positions(query: str, reference: str, ref_positions: List[int]) -> str:
    """ Extracts the given positions from a query alignment. The positions are
        adjusted to account for any gaps in the reference sequence.

        Arguments:
            query: the aligned query
            reference: the aligned reference
            ref_positions: the positions of interest in the unaligned reference

        Returns:
            a string containing the sequence elements at the adjusted reference
            positions
    """
    # adjust position of interest to account for gaps in the ref sequence alignment
    positions = []
    position_skipping_gaps = 0
    for i, amino in enumerate(reference):
        if amino == "-":
            continue
        if position_skipping_gaps in ref_positions:
            positions.append(i)
        position_skipping_gaps += 1
    # extract positions from query sequence
    return "".join([query[i] for i in positions])
