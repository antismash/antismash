# License: GNU Affero General Public License v3 or later
# A copy of GNU AGPL v3 should have been included in this software package in
# LICENSE.txt.

""" DNA encoding and promoter motif counting for fungal boundaries """

from dataclasses import dataclass
from enum import StrEnum
from functools import cache
from itertools import product

from Bio.Seq import reverse_complement
import numpy as np

from . import config as config


class MotifType(StrEnum):
    """ The DNA arrangement represented by a motif """

    INVERTED = "inverted"
    DIRECT = "direct"
    LINEAR = "linear"


@dataclass
class MotifCounts:
    """ Sorted unique motif IDs and their counts in one promoter """

    motif_ids: np.ndarray
    counts: np.ndarray

    @classmethod
    def from_counts(cls, counts: np.ndarray, id_offset: int = 0,
                    capped_max_count: int | None = 3) -> "MotifCounts":
        """ Selects nonzero counts, optionally caps them and assigns motif IDs

            Arguments:
                counts: the raw counts for this motif type, indexed by motif ID
                        within the type's own block
                id_offset: the start of this motif type's block in the shared
                           ID range
                capped_max_count: the maximum count kept for each motif in this
                                  promoter, or None to keep uncapped counts

            Returns:
                a MotifCounts instance with only the IDs whose counts are
                nonzero, optionally capped at capped_max_count
        """
        motif_ids = np.flatnonzero(counts)
        nonzero_counts = counts[motif_ids]
        if capped_max_count is not None:
            nonzero_counts = np.minimum(nonzero_counts, capped_max_count)
        return cls(
            motif_ids=motif_ids + id_offset,
            counts=nonzero_counts,
        )


def get_motif_type(motif_id: int) -> MotifType:
    """ Returns the repeat type that the given motif ID belongs to

        Each motif type occupies its own block of the shared ID range:
        inverted repeats start at 0, direct repeats start at
        config.DIRECT_MOTIF_ID_START and linear motifs start at
        config.LINEAR_MOTIF_ID_START.

        Arguments:
            motif_id: the ID of the motif to classify

        Returns:
            MotifType.INVERTED, MotifType.DIRECT or MotifType.LINEAR
    """
    if motif_id < config.DIRECT_MOTIF_ID_START:
        return MotifType.INVERTED
    if motif_id < config.LINEAR_MOTIF_ID_START:
        return MotifType.DIRECT
    return MotifType.LINEAR


def code_to_word(code: int, length: int) -> str:
    """ Converts a numerical k-mer code back into the DNA word it stands for.

        The code is read as a base-4 number with the digits mapped to bases in
        alphabetical order (A=0, C=1, G=2, T=3), inverting the IDs built by
        _kmer_ids. Codes needing fewer digits than the given length are padded
        with leading A's, so code 6 with length 3 gives "ACG".

        Arguments:
            code: the k-mer ID to convert, from 0 to 4 ** length - 1
            length: the number of bases in the resulting word, config.ARM for
                    repeat arms or config.LINEAR_LENGTH for linear motifs

        Returns:
            a string of that many bases over A, C, G and T
    """
    code_base4 = np.base_repr(code, base=4).zfill(length)
    return code_base4.translate(str.maketrans("0123", "ACGT"))


def motif_int_to_text(motif_id: int) -> str:
    """ Converts a motif ID into the sequence pattern that it represents

        Which pattern is built depends on the block the ID falls in, as
        reported by get_motif_type. Repeat patterns name both arms and the
        number of spacer bases between them, while a linear motif is a
        single word.

        Arguments:
            motif_id: the ID of the motif to convert

        Returns:
            a string of the form 'ACG-N5-CGT' for an inverted repeat, whose
            second arm is the reverse complement of the first, 'ACG-N5-ACG'
            for a direct repeat, whose arms are identical, or a plain word
            such as 'ACGTAC' for a linear motif
    """
    if motif_id >= config.LINEAR_MOTIF_ID_START:
        return code_to_word(motif_id - config.LINEAR_MOTIF_ID_START,
                            config.LINEAR_LENGTH)
    if motif_id >= config.DIRECT_MOTIF_ID_START:
        code, spacer = divmod(motif_id - config.DIRECT_MOTIF_ID_START,
                              config.SPACERS)
        arm = code_to_word(code, config.ARM)
        return f"{arm}-N{spacer}-{arm}"
    code, spacer = divmod(motif_id, config.SPACERS)
    arm = code_to_word(code, config.ARM)
    return f"{arm}-N{spacer}-{reverse_complement(arm)}"


@cache
def _kmer_ids(kmer_size: int) -> tuple[dict[str, int], np.ndarray]:
    """ Builds the k-mer IDs and their reverse-complement codes together

        IDs start at 0 in alphabetical order, e.g. AAA=0 and AAC=1 for size 3.
        Both lookups are cached per k-mer size and so every caller gets the
        same two objects and neither may be modified.

        Arguments:
            kmer_size: the number of bases in each k-mer

        Returns:
            a dictionary mapping k-mer strings to IDs, and an int64 array
            mapping each ID to its reverse-complement ID. For size 3 the
            dictionary maps "ACG" to 6 and the array maps 6 to 27 ("CGT")
    """
    kmer_ids = {"".join(kmer): number
                for number, kmer in enumerate(product("ACGT",
                                                      repeat=kmer_size))}
    reverse_codes = np.array(
        [kmer_ids[reverse_complement(kmer)] for kmer in kmer_ids],
        dtype=np.int64,
    )
    return kmer_ids, reverse_codes


def encode_kmers(promoter_seq: str,
                 kmer_size: int
                 ) -> tuple[np.ndarray, np.ndarray, np.ndarray]:
    """ Encodes overlapping k-mers, flags valid bases and encodes their
        reverse complements

        Every distinct k-mer has one fixed ID shared across all promoters.

        Arguments:
            promoter_seq: the promoter DNA sequence to encode
            kmer_size: the length of each k-mer, config.ARM for repeat arms
                       or config.LINEAR_LENGTH for linear motifs

        Returns:
            kmer_encodings: an int64 array with one ID per overlapping k-mer,
                            where entry i is the k-mer starting at base i
            valid_kmers: a bool array of the same length, True where the
                         k-mer contains only A, C, G or T. Entries where
                         valid_kmers is False still carry an ID(with every
                         non-ACGT base replaced by T) so callers must apply
                         the mask rather than relying on the IDs alone.
            numerical_form_of_reverse_complements:
                an int64 array of the IDs of the reverse-complement k-mers

        Example:
            For promoter_seq="ATGATG" and kmer_size=3:
                kmers = ["ATG", "TGA", "GAT", "ATG"]
                kmer_encodings = [14, 56, 35, 14]
                valid_kmers = [True, True, True, True]
                reverse_complements = ["CAT", "TCA", "ATC", "CAT"]
                numerical_form_of_reverse_complements = [19, 52, 13, 19]

            Each reverse complement stays at the same index as its original
            k-mer.
    """
    num_kmers = max(0, len(promoter_seq) - kmer_size + 1)
    kmer_encodings = np.zeros(num_kmers, dtype=np.int64)
    valid_kmers = np.ones(num_kmers, dtype=bool)
    if not num_kmers:
        return kmer_encodings, valid_kmers, np.zeros(0, dtype=np.int64)

    # create dict - key: kmer_seq, value: kmer_id
    kmer_to_ids_dict, reverse_codes = _kmer_ids(kmer_size)

    # loop over the promoter sequence
    for start in range(num_kmers):

        current_kmer_seq = promoter_seq[start:start + kmer_size]
        current_kmer_id = kmer_to_ids_dict.get(current_kmer_seq)

        if current_kmer_id is None:
            valid_kmers[start] = False
            # Unknown bases use T only to produce a placeholder ID.
            placeholder = "".join(base if base in "ACGT" else "T"
                                  for base in current_kmer_seq)
            current_kmer_id = kmer_to_ids_dict[placeholder]

        kmer_encodings[start] = current_kmer_id
    numerical_form_of_reverse_complements = reverse_codes[kmer_encodings]
    return kmer_encodings, valid_kmers, numerical_form_of_reverse_complements


def repeat_motif_counts(promoter_seq: str,
                        capped_max_count: int | None = config.MAX_COUNTS
                        ) -> MotifCounts:
    """ Counts the inverted and direct repeat motifs in the given sequence

        A repeat is two arms of config.ARM bases separated by a spacer of 0 to
        config.MAX_SPACER bases. In an inverted repeat the second arm is the
        reverse complement of the first; in a direct repeat the two arms are
        identical and get canonicalised. The spacer bases are never examined,
        so ACGNNCGT counts as an inverted repeat. Occurrences with an unknown
        base in either arm are skipped.

        Arguments:
            promoter_seq: the promoter sequence to scan, read 5' to 3'
            capped_max_count: the maximum count kept for each motif in this
                              promoter, or None to keep uncapped counts

        Returns:
            a MotifCounts instance with the sorted unique motif IDs and their
            counts, optionally capped at capped_max_count
    """

    # Skip promoters shorter than two adjacent repeat arms.
    if len(promoter_seq) < 2 * config.ARM:
        return MotifCounts.from_counts(np.zeros(0, dtype=np.int64),
                                       capped_max_count=capped_max_count)

    # see encode_kmers docstring for example of data
    (numerical_form_of_kmers, valid_kmers,
     numerical_form_of_reverse_complements) = encode_kmers(
         promoter_seq, config.ARM)

    # Each table stores counts at [arm_id, gap_size].
    num_arm_ids = 4 ** config.ARM
    inverted_counts = np.zeros((num_arm_ids, config.SPACERS), dtype=np.int64)
    direct_counts = np.zeros_like(inverted_counts)

    # Check every possible arm pair for each allowed gap size.
    for gap_size in range(config.MAX_SPACER + 1):
        second_arm_start_pos = config.ARM + gap_size
        first_arm_ids = numerical_form_of_kmers[:-second_arm_start_pos]

        # Keep only positions where both arms have known bases
        usable = (valid_kmers[:-second_arm_start_pos]
                  & valid_kmers[second_arm_start_pos:])

        # Find inverted repeats
        is_inverted = first_arm_ids == numerical_form_of_reverse_complements[
            second_arm_start_pos:]
        inverted_arm_ids = first_arm_ids[usable & is_inverted]

        # Find direct repeats (both arms are identical)
        is_direct = usable & (
            first_arm_ids == numerical_form_of_kmers[second_arm_start_pos:])
        direct_arm_ids = np.minimum(
            first_arm_ids[is_direct],
            numerical_form_of_reverse_complements[:-second_arm_start_pos][is_direct],
        )
        # Count occurrences of each arm for this gap size.
        inverted_counts[:, gap_size] = np.bincount(inverted_arm_ids,
                                                   minlength=num_arm_ids)
        direct_counts[:, gap_size] = np.bincount(direct_arm_ids,
                                                 minlength=num_arm_ids)

    # Flatten each arm's gap counts in order, with inverted repeats before
    # direct repeats. This preserves the motif IDs used by the scoring pipeline
    return MotifCounts.from_counts(
        np.concatenate((inverted_counts.ravel(), direct_counts.ravel())),
        capped_max_count=capped_max_count,
    )


def linear_motif_counts(promoter_seq: str,
                        capped_max_count: int | None = config.MAX_COUNTS
                        ) -> MotifCounts:
    """ Counts the linear motifs of config.LINEAR_LENGTH bases in the sequence

        A motif and its reverse complement share one ID, keeping the lower one.
        Occurrences with an unknown base are skipped.

        Arguments:
            promoter_seq: the promoter sequence to scan
            capped_max_count: the maximum count kept for each motif in this
                              promoter, or None to keep uncapped counts

        Returns:
            a MotifCounts instance with the sorted unique motif IDs and their
            counts, optionally capped at capped_max_count
    """

    if len(promoter_seq) < config.LINEAR_LENGTH:
        return MotifCounts.from_counts(np.zeros(0, dtype=np.int64),
                                       config.LINEAR_MOTIF_ID_START,
                                       capped_max_count=capped_max_count)
    (numerical_form_of_kmers, valid_kmers,
     numerical_form_of_reverse_complements) = encode_kmers(
         promoter_seq, config.LINEAR_LENGTH)
    # A linear motif and its reverse complement are counted as one motif.
    canonical_ids = np.minimum(numerical_form_of_kmers,
                               numerical_form_of_reverse_complements)
    counts = np.bincount(canonical_ids[valid_kmers])
    return MotifCounts.from_counts(counts, config.LINEAR_MOTIF_ID_START,
                                   capped_max_count=capped_max_count)
