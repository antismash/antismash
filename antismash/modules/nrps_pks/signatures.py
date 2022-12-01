# License: GNU Affero General Public License v3 or later
# A copy of GNU AGPL v3 should have been included in this software package in LICENSE.txt.

"""
Handle extracting the 10 AA and 34 AA signatures of NRPS adenylation domains.
"""

from antismash.common import path, subprocessing
from antismash.detection.nrps_pks_domains import ModularDomain

REF_SEQUENCE = "P0C062_A1"
A34_POSITIONS_FILENAME = path.get_full_path(__file__, "external", "NRPSPredictor2", "A34positions.txt")
APOSITION_FILENAME = path.get_full_path(__file__, "external", "NRPSPredictor2", "Apositions.txt")
ILLEGAL_CHARS = "!@#$%^&*(){}:\"<>?/.,';][`~1234567890*-+-=_\\|"
ADOMAINS_FILENAME = path.get_full_path(__file__, "external", "NRPSPredictor2", "A_domains_muscle.fasta")
START_POSITION = 66


def read_positions(filename: str, start_position: int) -> list[int]:
    """ Loads positions from a tab-separated file. Positions are relative to the start_position.

        Arguments:
            filename: the path to the file containing the positions
            start_position: a relative start position to adjust all positions by

        Returns:
            a list of ints, one for each position found in the file
    """
    with open(filename, "r", encoding="utf-8") as data:
        text = data.read().strip()
    results = []
    for i in text.split("\t"):
        results.append(int(i) - start_position)
    return results


def build_position_list(positions: list[int], reference_seq: str) -> list[int]:
    """ Adjusts a list of positions to account for gaps in the reference sequence

        Arguments:
            positions: a list of ints that represent positions of interest in
                       the reference sequence
            reference_seq: the (aligned) reference sequence

        Returns:
            a new list of positions, each >= the original position
    """
    poslist = []
    position = 0
    for i, ref in enumerate(reference_seq):
        if ref != "-":
            if position in positions:
                poslist.append(i)
            position += 1
    return poslist


def verify_good_sequence(sequence: str) -> bool:
    """ Ensures a sequence is valid """
    for char in ILLEGAL_CHARS:
        if char in sequence:
            return False
    return True


def extract(sequence: str, positions: list[int]) -> str:
    """ Extracts a signature from an aligned sequence based on the provided
        positions. Accounts for gaps by looking behind or, if behind is already
        in the position list, ahead.

        Arguments:
            sequence: the aligned sequence to extract a signature from
            positions: the list of positions within the sequence to use

        Returns:
            the extracted signature as a string
    """
    seq = []
    for position in positions:
        aa = sequence[position]
        if aa == "-":
            if position - 1 not in positions:
                aa = sequence[position - 1]
            elif position + 1 not in positions:
                aa = sequence[position + 1]
        seq.append(aa)
    return "".join(seq)


def get_a_dom_signatures(domain: ModularDomain) -> tuple[str, str]:
    """ Extract 10 / 34 AA NRPS signatures from A domains """
    assert " " not in domain.get_name()
    assert verify_good_sequence(domain.translation)
    # Run muscle and collect sequence positions from file
    alignments = subprocessing.run_muscle_single(domain.get_name(), domain.translation, ADOMAINS_FILENAME)

    domain_alignment = alignments[domain.get_name()]
    reference_alignment = alignments[REF_SEQUENCE]

    positions = read_positions(APOSITION_FILENAME, START_POSITION)
    # Count residues in ref sequence and put positions in list
    poslist = build_position_list(positions, reference_alignment)

    # Extract positions from query sequence
    query_sig_seq = extract(domain_alignment, poslist)
    # Add fixed lysine 517
    query_sig_seq += "K"

    # repeat with 34 AA codes
    angpositions = read_positions(A34_POSITIONS_FILENAME, START_POSITION)
    poslist = build_position_list(angpositions, reference_alignment)

    return query_sig_seq, extract(domain_alignment, poslist)
