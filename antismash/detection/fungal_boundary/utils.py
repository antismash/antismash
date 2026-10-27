# License: GNU Affero General Public License v3 or later
# A copy of GNU AGPL v3 should have been included in this software package in
# LICENSE.txt.

"""DNA encoding, motif detection and window scoring for fungal boundaries."""

from dataclasses import dataclass
from functools import cache
from itertools import product

from Bio.Seq import reverse_complement
import numpy as np

from . import config


@dataclass
class MotifCounts:
    """Sorted unique motif IDs and their capped counts in one promoter."""

    motif_ids: np.ndarray
    counts: np.ndarray

    @classmethod
    def from_counts(cls, counts: np.ndarray, id_offset: int = 0) -> "MotifCounts":
        """Keep nonzero counts, cap them and convert indices to motif IDs.

        id_offset locates this motif type's block in the shared ID range.
        """
        motif_ids = np.flatnonzero(counts)
        return cls(
            motif_ids=motif_ids + id_offset,
            counts=np.minimum(counts[motif_ids], config.MAX_COUNTS),
        )


@dataclass
class WindowScores:
    """Each motif's best score, hit counts and inclusive window endpoints.

    Each array is indexed by motif ID. observed holds the window's total
    hits, and with_site counts its genes carrying at least one hit. left and
    right are inclusive row indices in the gene-by-motif count table.
    """

    llr: np.ndarray
    observed_hits: np.ndarray
    genes_with_motif: np.ndarray
    left: np.ndarray
    right: np.ndarray


def running_totals(values: np.ndarray) -> np.ndarray:
    """Sum each row from window start to window end, with an initial
    zero column.


    Arguments:
        values: one row per motif, one column per gene; a boolean array is
            summed as 0 and 1, which is how genes carrying a motif are counted

    Returns:
        A float array with one more column than values, where entry [i, j] is
        the sum of the first j values of row i.
    """
    totals = np.zeros((len(values), values.shape[1] + 1))
    values.cumsum(axis=1, out=totals[:, 1:])
    return totals


def poisson_scores(observed: np.ndarray, expected: np.ndarray,
                   eligible: np.ndarray) -> np.ndarray:
    """Calculate the Poisson score only for eligible a window.

        Arguments:
        observed: hits counted in each window
        expected: hits the background rate predicts for each window; must be
            strictly positive wherever eligible is True
        eligible: which windows to score, already checked by the caller for
            enrichment and support

    Returns:
        An array shaped like observed, holding each eligible window's score
        and 0 everywhere else.
    """
    scores = np.zeros(observed.shape)
    hits = observed[eligible]
    background_rate = expected[eligible]
    scores[eligible] = hits * np.log(hits / background_rate) - (hits - background_rate)
    return scores


def _initial_window_scores(counts: np.ndarray, start: int, end: int) -> WindowScores:
    """Default every motif to a zero score and the first window supplied."""
    first_window = counts[start:end + 1]
    n_motifs = counts.shape[1]
    return WindowScores(
        llr=np.zeros(n_motifs),
        observed_hits=first_window.sum(axis=0),
        genes_with_motif=(first_window > 0).sum(axis=0, dtype=float),
        left=np.full(n_motifs, start, dtype=int),
        right=np.full(n_motifs, end, dtype=int),
    )


def _update_best_windows(best: WindowScores, motif_ids: np.ndarray,
                         start: int, ends: np.ndarray, window_llr: np.ndarray,
                         observed_hits: np.ndarray, genes_with_motif: np.ndarray) -> None:
    """Keep each motif's best window from this batch when its score improves."""
    # argmax chooses the first window end when scores are tied.
    best_column = window_llr.argmax(axis=1)
    batch_rows = np.arange(len(motif_ids))
    best_scores = window_llr[batch_rows, best_column]
    better = best_scores > best.llr[motif_ids]
    winners = motif_ids[better]
    winning_rows = batch_rows[better]
    winning_columns = best_column[better]
    best.llr[winners] = best_scores[better]
    best.observed_hits[winners] = observed_hits[winning_rows, winning_columns]
    best.genes_with_motif[winners] = genes_with_motif[winning_rows, winning_columns]
    best.left[winners] = start
    best.right[winners] = ends[winning_columns]


def find_best_windows(counts: np.ndarray, rate: np.ndarray,
                      window_starts: np.ndarray, window_ends: np.ndarray,
                      center: int) -> WindowScores:
    """Find the highest-scoring window for every motif.

    A window scores only when it holds at least config.MIN_GENES_WITH_MOTIFS
    genes carrying the motif and more hits than the background rate predicts;
    the score itself comes from poisson_scores. Anything else keeps a score
    of 0.

    Windows are grouped by their start window, so a single pass covers every end
    for that start, for every motif at once. running_totals reduces each window's
    hit count to one subtraction, and motifs present in too few genes to ever
    qualify are dropped before any scoring.

    A window whose first or last gene lacks the motif can never be the best
    one: dropping that gene leaves the same hits over fewer genes, which
    scores higher, and the shorter window is supplied separately. Such
    windows are skipped, unless the empty gene is the anchor, which cannot be
    dropped, or the window is already at config.MIN_GENES.

    Arguments:
        counts: motif counts per gene, genes down the rows and motif IDs
            across the columns
        rate: each motif's background hits per gene, indexed by motif ID
        window_starts: each window's first row in counts
        window_ends: each window's last row, inclusive. Both arrays must be
            ordered by start then end, with consecutive ends for each start
        center: the anchor's row in counts

    Returns:
        One entry per column of counts: the best window's score, its hits,
        how many of its genes carry the motif, and its inclusive endpoints as
        rows of counts. Tied scores keep the window met first, which is the
        smallest start and then the smallest end. A motif that never scored
        keeps 0 and the first window supplied, so check llr before reading a
        motif's endpoints.
    """
    # Step 1: start every motif with no score and the first window supplied,
    # which is what a motif that never scores keeps.
    best = _initial_window_scores(counts, window_starts[0], window_ends[0])

    # Step 2: keep only the genes that can belong to a window, put one motif
    # per row, and take running totals so that any window's hits and its
    # supporting genes each cost one subtraction. Motifs occurring in too few
    # genes to ever qualify are dropped here and keep the step 1 defaults.
    first_row = window_starts[0]
    end_row = window_ends.max() + 1
    # Transpose the matrix so that it is easier to work with motifs
    counts_by_motif = counts[first_row:end_row].T
    has_motif = counts_by_motif > 0
    scorable_motifs = np.flatnonzero(has_motif.sum(axis=1) >= config.MIN_GENES_WITH_MOTIFS)
    if scorable_motifs.size == 0:
        return best
    counts_by_motif = counts_by_motif[scorable_motifs]
    has_motif = has_motif[scorable_motifs]
    hit_totals = running_totals(counts_by_motif)
    with_motif_totals = running_totals(has_motif)

    # Step 3: score every window beginning at this start, for all the motifs
    # still in play, in one pass.
    def score_start(start: int, ends: np.ndarray) -> None:
        """Score all window ends for one start using the shared motif totals."""
        window_sizes = ends - start + 1
        start_column = start - first_row
        first_end_column = ends[0] - first_row
        stop_end_column = ends[-1] - first_row + 1

        # A minimum-size window may need an empty end gene to stay large enough.
        at_minimum_size = (window_sizes[0] == config.MIN_GENES
                           and config.MIN_GENES >= config.MIN_GENES_WITH_MOTIFS)
        if start == center or at_minimum_size:
            scored_rows: np.ndarray = np.arange(len(scorable_motifs))
        else:
            scored_rows = np.flatnonzero(has_motif[:, start_column])
        if scored_rows.size == 0:
            return

        # Rows are motifs, columns are window ends. None adds a column axis so
        # each motif's starting total and rate apply to every window end.
        observed_hits = (
            hit_totals[scored_rows, first_end_column + 1:stop_end_column + 1]
            - hit_totals[scored_rows, start_column, None])
        genes_with_motif = (
            with_motif_totals[scored_rows, first_end_column + 1:stop_end_column + 1]
            - with_motif_totals[scored_rows, start_column, None])
        expected_hits = rate[scorable_motifs[scored_rows], None] * window_sizes

        start_supported = (has_motif[scored_rows, start_column, None]
                           | (start == center) | (window_sizes == config.MIN_GENES))
        end_supported = (has_motif[scored_rows, first_end_column:stop_end_column]
                         | (ends == center) | (window_sizes == config.MIN_GENES))
        eligible = ((genes_with_motif >= config.MIN_GENES_WITH_MOTIFS)
                    & (observed_hits > expected_hits)
                    & start_supported & end_supported)
        window_llr = poisson_scores(observed_hits, expected_hits, eligible)

        # Step 4: take each motif's best window from this pass and keep it
        # where it beats what an earlier start found.
        _update_best_windows(best, scorable_motifs[scored_rows], start, ends,
                             window_llr, observed_hits, genes_with_motif)

    for start in np.unique(window_starts):
        score_start(start, window_ends[window_starts == start])
    return best


def code_to_word(code: int, length: int) -> str:
    """Spell out the DNA word that a single base-4 code stands for.

    The inverse of word_to_code, and of one entry of encode_kmers: the code
    is read as a base-4 number with A=0, C=1, G=2 and T=3 and padded to
    length bases.

    Arguments:
        code: the word's base-4 value, from 0 to 4 ** length - 1
        length: how many bases to spell out, config.ARM for a repeat arm or
                config.LINEAR_LENGTH for a linear motif

    Returns:
        The word as uppercase DNA text, always length bases long.
    """
    to_dna = str.maketrans("0123", "ACGT")
    code_base4 = np.base_repr(code, base=4)
    code_padded = code_base4.zfill(length)
    trimmed = code_padded[-length:]
    dna = trimmed.translate(to_dna)
    return dna


def word_to_code(sequence: str) -> int:
    """Read a whole DNA word as a single base-4 number.

    The inverse of code_to_word, and the scalar form of encode_kmers: the
    word becomes one number, so "ACG" gives 0 * 16 + 1 * 4 + 2 = 6.
    Only uppercase A, C, G and T are accepted; any other character,
    and the empty string, raise ValueError.

    Arguments:
        sequence: a DNA word, A, C, G and T only

    Returns:
        The word's base-4 value, from 0 to 4 ** len(sequence) - 1.
    """
    to_digits = str.maketrans("ACGT", "0123")
    code = int(sequence.translate(to_digits), 4)
    return code


def motif_type(motif_id: int) -> str:
    """Name the kind of motif that an ID belongs to.

    The ID alone decides, because each kind owns its own block of the ID
    range; see config.MOTIF_TYPE_RANGE.

    Arguments:
        motif_id: a motif ID, from 0 to config.N_MOTIFS - 1

    Returns:
        "inverted", "direct" or "linear".
    """
    if motif_id >= config.LINEAR_START:
        return "linear"
    if motif_id >= config.DIRECT_START:
        return "direct"
    return "inverted"


def motif_int_to_text(motif_id: int) -> str:
    """Converts a motif id into a .

    Arguments:
        motif_id: a motif ID, from 0 to config.N_MOTIFS - 1

    Returns:
        The motif as DNA text.
    """
    if motif_id >= config.LINEAR_START:
        return code_to_word(motif_id - config.LINEAR_START, config.LINEAR_LENGTH)
    if motif_id >= config.DIRECT_START:
        code, spacer = divmod(motif_id - config.DIRECT_START, config.SPACERS)
        arm = code_to_word(code, config.ARM)
        return f"{arm}-N{spacer}-{arm}"
    code, spacer = divmod(motif_id, config.SPACERS)
    arm = code_to_word(code, config.ARM)
    return f"{arm}-N{spacer}-{reverse_complement(arm)}"


@cache
def _kmer_ids(kmer_size: int) -> tuple[dict[str, int], np.ndarray]:
    """Build k-mer IDs and their reverse-complement codes together.

    IDs start at 0 in alphabetical order (e.g. AA=0 and AC=1 for size 2).
    Both lookups are cached per k-mer size and reused across promoters.

    Arguments:
        kmer_size: number of bases in each k-mer

    Returns:
        A dictionary mapping k-mer strings to IDs, and an int64 array
        mapping each ID to its reverse-complement ID. For size 3,
        the dictionary maps "ACG" to 6 and the array maps 6 to 27 ("CGT").
    """
    kmer_ids = {"".join(kmer): number
                for number, kmer in enumerate(product("ACGT", repeat=kmer_size))}
    reverse_codes = np.array(
        [kmer_ids[reverse_complement(kmer)] for kmer in kmer_ids],
        dtype=np.int64,
    )
    return kmer_ids, reverse_codes


def encode_kmers(promotor_seq: str,
                 kmer_size: int) -> tuple[np.ndarray, np.ndarray, np.ndarray]:
    """Encode overlapping k-mers, flag valid bases and encode reverse complements.

    Every distinct k-mer has one fixed ID shared across all promoters.

    Arguments:
        promotor_seq: promoter DNA sequence
        kmer_size: length of each k-mer, config.ARM for repeat arms or
                   config.LINEAR_LENGTH for linear motifs
    Returns:
        kmer_encodings: int64 array, one number per overlapping k-mer, where
                        entry i is the k-mer starting at base i
        valid_kmers: bool array of the same length, True where the k-mer
                     contains only A, C, G or T
        numerical_form_of_reverse_complements:
            int64 array of IDs which correspond to the reverse-complement kmers.

    Example:
        For promotor_seq="ATGATG" and kmer_size=3:
            kmers = ["ATG", "TGA", "GAT", "ATG"]
            numerical_form_of_kmers = [14, 56, 35, 14]
            valid_kmers = [True, True, True, True]
            reverse_complements = ["CAT", "TCA", "ATC", "CAT"]
            numerical_form_of_reverse_complements = [19, 52, 13, 19]

        Each reverse complement stays at the same index as its original k-mer.
    """
    num_kmers = max(0, len(promotor_seq) - kmer_size + 1)
    kmer_encodings = np.zeros(num_kmers, dtype=np.int64)
    valid_kmers = np.ones(num_kmers, dtype=bool)
    if not num_kmers:
        return kmer_encodings, valid_kmers, np.zeros(0, dtype=np.int64)

    # create dict - key: kmer_seq, value: kmer_id
    kmer_to_ids_dict, reverse_codes = _kmer_ids(kmer_size)

    # loop over the promoter sequence
    for start in range(num_kmers):

        current_kmer_seq = promotor_seq[start:start + kmer_size]
        current_kmer_id = kmer_to_ids_dict.get(current_kmer_seq)

        if current_kmer_id is None:
            valid_kmers[start] = False
            # Unknown bases use T only to produce a placeholder ID.
            placeholder = "".join(base if base in "ACGT" else "T" for base in current_kmer_seq)
            current_kmer_id = kmer_to_ids_dict[placeholder]

        kmer_encodings[start] = current_kmer_id
    numerical_form_of_reverse_complements = reverse_codes[kmer_encodings]
    return kmer_encodings, valid_kmers, numerical_form_of_reverse_complements


def repeat_motif_counts(promotor_seq: str) -> MotifCounts:
    """Count inverted and direct repeat motifs in the sequence.

    A repeat is two arms of config.ARM bases separated by a spacer of 0 to
    config.MAX_SPACER bases. In an inverted repeat the second arm is the
    reverse complement of the first; in a direct repeat the two arms are
    identical. The spacer bases are never examined, so ACGNNCGT counts as
    an inverted repeat. Occurrences with an unknown base in either arm
    are skipped.

    Arguments:
            promotor_seq: the promoter sequence to scan, read 5' to 3'

    Returns:
            Sorted unique motif IDs and counts capped at config.MAX_COUNTS.
    """

    # Skip promoters shorter than two adjacent repeat arms.
    if len(promotor_seq) < 2 * config.ARM:
        return MotifCounts.from_counts(np.zeros(0, dtype=np.int64))

    # see encode_kmers docstring for example of data
    numerical_form_of_kmers, valid_kmers, numerical_form_of_reverse_complements = encode_kmers(
        promotor_seq, config.ARM)

    # Each table stores counts at [arm_id, gap_size].
    num_arm_ids = 4 ** config.ARM
    inverted_counts = np.zeros((num_arm_ids, config.SPACERS), dtype=np.int64)
    direct_counts = np.zeros_like(inverted_counts)

    # Check every possible arm pair for each allowed gap size.
    for gap_size in range(config.MAX_SPACER + 1):
        second_arm_start_pos = config.ARM + gap_size
        first_arm_ids = numerical_form_of_kmers[:-second_arm_start_pos]

        # Keep only positions where both arms have known bases
        usable = valid_kmers[:-second_arm_start_pos] & valid_kmers[second_arm_start_pos:]

        # Find inverted repeats
        is_inverted = first_arm_ids == numerical_form_of_reverse_complements[
            second_arm_start_pos:]
        inverted_arm_ids = first_arm_ids[usable & is_inverted]

        # Find direct repeats (both arms are identical)
        is_direct = usable & (first_arm_ids == numerical_form_of_kmers[second_arm_start_pos:])
        direct_arm_ids = np.minimum(
            first_arm_ids[is_direct],
            numerical_form_of_reverse_complements[:-second_arm_start_pos][is_direct],
        )

        # Count occurrences of each arm for this gap size.
        inverted_counts[:, gap_size] = np.bincount(inverted_arm_ids, minlength=num_arm_ids)
        direct_counts[:, gap_size] = np.bincount(direct_arm_ids, minlength=num_arm_ids)

    # Flatten each arm's gap counts in order, with inverted repeats before
    # direct repeats. This preserves the motif IDs used by the scoring pipeline.
    return MotifCounts.from_counts(np.concatenate((inverted_counts.ravel(), direct_counts.ravel())))


def linear_motif_counts(promotor_seq: str) -> MotifCounts:
    """Count linear motifs of config.LINEAR_LENGTH bases in the sequence.

    A motif and its reverse complement share one ID. Occurrences with an
    unknown base are skipped.

        Arguments:
        promotor_seq: the promoter sequence to scan

    Returns:
        Sorted unique motif IDs and counts capped at config.MAX_COUNTS.
    """

    if len(promotor_seq) < config.LINEAR_LENGTH:
        return MotifCounts.from_counts(np.zeros(0, dtype=np.int64), config.LINEAR_START)
    numerical_form_of_kmers, valid_kmers, numerical_form_of_reverse_complements = encode_kmers(
        promotor_seq, config.LINEAR_LENGTH)
    # A linear motif and its reverse complement are counted as one motif.
    canonical_ids = np.minimum(numerical_form_of_kmers, numerical_form_of_reverse_complements)
    counts = np.bincount(canonical_ids[valid_kmers])
    return MotifCounts.from_counts(counts, config.LINEAR_START)
