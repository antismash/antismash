# License: GNU Affero General Public License v3 or later
# A copy of GNU AGPL v3 should have been included in this software package in LICENSE.txt.

"""DNA encoding, motif detection and window scoring for fungal boundaries."""

from dataclasses import dataclass

import numpy as np

from .config import (
    ARM, DIRECT_START, LINEAR_LENGTH, LINEAR_START, MIN_GENES,
    MIN_GENES_WITH_MOTIFS, SPACERS,
)


@dataclass
class WindowScores:
    """Each motif's best score, hit counts and inclusive window endpoints."""

    llr: np.ndarray
    observed: np.ndarray
    with_site: np.ndarray
    left: np.ndarray
    right: np.ndarray


def running_totals(values: np.ndarray) -> np.ndarray:
    """Sum each row from left to right, with an initial zero column.

    A window's sum is then totals[:, right + 1] - totals[:, left].
    """
    totals = np.zeros((len(values), values.shape[1] + 1))
    values.cumsum(axis=1, out=totals[:, 1:])
    return totals


def poisson_scores(observed: np.ndarray, expected: np.ndarray,
                   eligible: np.ndarray) -> np.ndarray:
    """Calculate the Poisson score only for eligible, enriched windows."""
    scores = np.zeros(observed.shape)
    hits = observed[eligible]
    background = expected[eligible]
    scores[eligible] = hits * np.log(hits / background) - (hits - background)
    return scores


def score_motif_windows(counts: np.ndarray, rate: np.ndarray,
                        lefts: np.ndarray, rights: np.ndarray,
                        center: int) -> WindowScores:
    """Find each motif's best window, processing one left endpoint at a time.

    Supply consecutive right endpoints for each left endpoint, ordered by
    left then right. Equal scores keep the first window, including when every
    score for a motif is zero.

    Empty end genes can be trimmed to increase enrichment, unless the anchor
    or MIN_GENES prevents trimming. Only those possible winners need scoring.
    """
    first_window = counts[lefts[0]:rights[0] + 1]
    n_motifs = counts.shape[1]
    best = WindowScores(
        llr=np.zeros(n_motifs),
        observed=first_window.sum(axis=0),
        with_site=(first_window > 0).sum(axis=0, dtype=float),
        left=np.full(n_motifs, lefts[0], dtype=int),
        right=np.full(n_motifs, rights[0], dtype=int),
    )

    # Keep only genes that can belong to a window, with genes along each row.
    first = lefts[0]
    last = rights.max() + 1
    motif_counts = counts[first:last].T.copy()
    present = motif_counts > 0
    motifs = np.flatnonzero(present.sum(axis=1) >= MIN_GENES_WITH_MOTIFS)
    if not len(motifs):
        return best
    motif_counts = motif_counts[motifs]
    present = present[motifs]
    hits = running_totals(motif_counts)
    sites = running_totals(present)

    for left in np.unique(lefts):
        right = rights[lefts == left]
        size = right - left + 1
        left_column = left - first
        start = right[0] - first
        stop = right[-1] - first + 1

        # A minimum-size window may need an empty end gene to stay large enough.
        minimum_window = (size[0] == MIN_GENES
                          and MIN_GENES >= MIN_GENES_WITH_MOTIFS)
        if left == center or minimum_window:
            rows = np.arange(len(motifs))
        else:
            rows = np.flatnonzero(present[:, left_column])
        if not len(rows):
            continue

        # Rows are motifs, columns are right endpoints. None adds a column axis
        # so each motif's starting total and rate apply to every right endpoint.
        observed = hits[rows, start + 1:stop + 1] - hits[rows, left_column, None]
        with_site = sites[rows, start + 1:stop + 1] - sites[rows, left_column, None]
        expected = rate[motifs[rows], None] * size

        left_supported = (present[rows, left_column, None]
                          | (left == center) | (size == MIN_GENES))
        right_supported = (present[rows, start:stop]
                           | (right == center) | (size == MIN_GENES))
        eligible = ((with_site >= MIN_GENES_WITH_MOTIFS) & (observed > expected)
                    & left_supported & right_supported)
        llr = poisson_scores(observed, expected, eligible)

        # argmax chooses the first right endpoint when scores are tied.
        column = llr.argmax(axis=1)
        row = np.arange(len(rows))
        scores = llr[row, column]
        better = scores > best.llr[motifs[rows]]
        winners = motifs[rows[better]]
        row = row[better]
        column = column[better]
        best.llr[winners] = scores[better]
        best.observed[winners] = observed[row, column]
        best.with_site[winners] = with_site[row, column]
        best.left[winners] = left
        best.right[winners] = right[column]
    return best


def reverse_complement(sequence: str) -> str:
    """ Returns reverse complement of the sequence"""
    complement = str.maketrans("ACGT", "TGCA")
    reverse_com = sequence.translate(complement)[::-1]
    return reverse_com


def num_to_dna(code: int, length: int) -> str:
    """converts numeric to DNA"""
    to_dna = str.maketrans("0123", "ACGT")
    code_base4 = np.base_repr(code, base=4)
    code_padded = code_base4.zfill(length)
    trimmed = code_padded[-length:]
    dna = trimmed.translate(to_dna)
    return dna


def dna_to_num(sequence: str) -> int:
    """converts DNA to numeric value"""
    to_digits = str.maketrans("ACGT", "0123")
    code = int(sequence.translate(to_digits), 4)
    return code


def reverse_num(length: int) -> np.ndarray:
    """Look up table for the conversion to a numeric value
    of a reverse complement"""
    codes = np.zeros(4 ** length, dtype=np.int64)
    for code in range(4 ** length):
        codes[code] = dna_to_num(reverse_complement(num_to_dna(code, length)))
    return codes


# Build lookup tables once at import time for reuse during motif detection.
_BASE_CODES = {ord(base): code for code, base in enumerate("ACGT")}
_BASE_TRANSLATION = bytes(_BASE_CODES.get(base, 4) for base in range(256))
REVERSE_ARM = reverse_num(ARM)
REVERSE_LINEAR = reverse_num(LINEAR_LENGTH)


def motif_type(motif_id: int) -> str:
    """Returns a type of motif using it's motif ID"""
    if motif_id >= LINEAR_START:
        return "linear"
    if motif_id >= DIRECT_START:
        return "direct"
    return "inverted"


def motif_int_to_text(motif_id: int) -> str:
    """Convert motif ID to a string motif"""
    if motif_id >= LINEAR_START:
        return num_to_dna(motif_id - LINEAR_START, LINEAR_LENGTH)
    if motif_id >= DIRECT_START:
        code, spacer = divmod(motif_id - DIRECT_START, SPACERS)
        arm = num_to_dna(code, ARM)
        return f"{arm}-N{spacer}-{arm}"
    code, spacer = divmod(motif_id, SPACERS)
    arm = num_to_dna(code, ARM)
    return f"{arm}-N{spacer}-{reverse_complement(arm)}"


def base_to_num(sequence: str) -> np.ndarray:
    """Encode A/C/G/T as 0/1/2/3 and every other character as 4."""
    # Convert the DNA text to bytes; replace unsupported characters with "?".
    encoded = sequence.encode("ascii", errors="replace")
    # Use the lookup table: A=0, C=1, G=2, T=3, and unknown bases=4.
    encoded = encoded.translate(_BASE_TRANSLATION)
    # Read the numbers as a NumPy array without copying them.
    bases = np.frombuffer(encoded, dtype=np.uint8)
    return bases


def kmer_to_num(bases: np.ndarray,
                kmer_size: int) -> tuple[np.ndarray, np.ndarray]:
    """Convert each overlapping k-mer to a number.

    Also return a boolean array: True means the k-mer contains only A, C, G or T.
    Process all k-mers together, one base position at a time.
    """
    num_kmers = max(0, len(bases) - kmer_size + 1)
    codes = np.zeros(num_kmers, dtype=np.int64)
    clean_kmers = np.ones(num_kmers, dtype=bool)
    base_codes = np.minimum(bases, 3)
    known_bases = bases != 4

    for offset in range(kmer_size):
        codes *= 4
        codes += base_codes[offset:offset + num_kmers]
        clean_kmers &= known_bases[offset:offset + num_kmers]
    return codes, clean_kmers


def repeat_motif_ids(sequence: str) -> np.ndarray:
    """Motif ID of every inverted and direct repeat in the sequence.
    A repeat is two arms of ARM bases separetad by a spacer. Each ocurrence
    givs one ID. Arms with uknown bases are skipped."""
    bases = base_to_num(sequence)
    if len(bases) < 2 * ARM:
        # checks that the promoter is at least the size of the smallest
        # possible repeat motif
        return np.zeros(0, dtype=np.int64)
    numerical_form_of_kmers, clean = kmer_to_num(bases, ARM)
    reverse_kmers = REVERSE_ARM[numerical_form_of_kmers]
    found = []
    for spacer in range(SPACERS):
        possible_motif = len(bases) - 2 * ARM - spacer + 1
        if possible_motif <= 0:
            break
        # Extract possible left arms
        left_arms = numerical_form_of_kmers[:possible_motif]
        # Extract possible right arms
        right_start = ARM + spacer
        right_arms = (numerical_form_of_kmers[right_start:right_start
                                              + possible_motif])

        # Keep onlu positions where both arms have known bases
        usable = (clean[:possible_motif]
                  & clean[right_start:right_start + possible_motif])

        # Find inverted repeats
        is_inverted = left_arms == reverse_kmers[right_start:right_start + possible_motif]
        inverted_codes = left_arms[usable & is_inverted]

        # Find direct repeats (both arms are identical)
        direct_codes = left_arms[usable & (left_arms == right_arms)]

        # Turn each repet into an ID that encodes the arm and spacer
        found.append(inverted_codes * SPACERS + spacer)
        found.append(DIRECT_START
                     + np.minimum(direct_codes, REVERSE_ARM[direct_codes])
                     * SPACERS + spacer)
    return np.concatenate(found)


def linear_motif_ids(sequence: str) -> np.ndarray:
    """Motif ID of every linera motif in the sequence.
    Each ocurrence gives one ID."""
    bases = base_to_num(sequence)
    if len(bases) < LINEAR_LENGTH:
        return np.zeros(0, dtype=np.int64)
    numerical_form_of_kmers, clean = kmer_to_num(bases, LINEAR_LENGTH)
    codes = numerical_form_of_kmers[clean]
    # A linear motif and its reverse complement are counted as one motif.
    return LINEAR_START + np.minimum(codes, REVERSE_LINEAR[codes])
