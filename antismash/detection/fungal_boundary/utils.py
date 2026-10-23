# License: GNU Affero General Public License v3 or later
# A copy of GNU AGPL v3 should have been included in this software package in
# LICENSE.txt.

""" DNA encoding and promoter motif counting for fungal boundaries """

from collections.abc import Mapping
from dataclasses import dataclass
from enum import StrEnum
from functools import cache
from itertools import product
from types import MappingProxyType

from Bio.Seq import reverse_complement
import numpy as np

from .config import config


class MotifType(StrEnum):
    """ The DNA arrangement represented by a motif """

    INVERTED = "inverted"
    DIRECT = "direct"
    LINEAR = "linear"


@dataclass
class MotifCounts:
    """ Unique motif IDs and their counts in one promoter """

    motif_ids: np.ndarray
    counts: np.ndarray

    @classmethod
    def from_counts(cls, counts: np.ndarray,
                    capped_max_count: int | None = config.MAX_MOTIFS_PER_PROMOTER
                    ) -> "MotifCounts":
        """ Selects motif IDs with nonzero counts and limits each count to
            capped_max_count.

            Arguments:
                counts: the raw counts, indexed by motif ID
                capped_max_count: the maximum count kept for each motif in this
                                  promoter, or None to keep uncapped counts

            Returns:
                a MotifCounts instance containing motif IDs and their counts
        """
        motif_ids = np.flatnonzero(counts)
        nonzero_counts = counts[motif_ids]
        if capped_max_count is not None:
            nonzero_counts = np.minimum(nonzero_counts, capped_max_count)
        return cls(
            motif_ids=motif_ids,
            counts=nonzero_counts,
        )


def get_motif_type(motif_id: int, *,
                   motif_type_id_boundaries: tuple[int, int] = config.motif_type_id_boundaries
                   ) -> MotifType:
    """ Returns the motif type that the given motif ID belongs to

        Arguments:
            motif_id: the ID of the motif to classify
            motif_type_id_boundaries: the first direct-repeat and linear motif IDs,
                                      in that order

        Returns:
            the motif type corresponding to the given motif ID
    """
    if motif_id < motif_type_id_boundaries[0]:
        return MotifType.INVERTED
    if motif_id < motif_type_id_boundaries[1]:
        return MotifType.DIRECT
    return MotifType.LINEAR


def _decode_kmer_id(code: int, length: int) -> str:
    """ Converts a k-mer ID into a DNA sequence of the given length """
    code_base4 = np.base_repr(code, base=4).zfill(length)
    return code_base4.translate(str.maketrans("0123", "ACGT"))


def convert_motif_id_to_text(motif_id: int) -> str:
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
    if motif_id >= config.linear_motif_id_start:
        return _decode_kmer_id(motif_id - config.linear_motif_id_start,
                               config.LINEAR_LENGTH)
    if motif_id >= config.direct_motif_id_start:
        code, spacer = divmod(motif_id - config.direct_motif_id_start,
                              config.spacers)
        arm = _decode_kmer_id(code, config.BASES_PER_ARM)
        return f"{arm}-N{spacer}-{arm}"
    code, spacer = divmod(motif_id, config.spacers)
    arm = _decode_kmer_id(code, config.BASES_PER_ARM)
    return f"{arm}-N{spacer}-{reverse_complement(arm)}"


@cache
def _get_kmer_to_ids(kmer_size: int) -> tuple[Mapping[str, int], np.ndarray]:
    """ Maps each k-mer sequence to an ID and each ID to its reverse-complement ID

        Arguments:
            kmer_size: the number of bases in each k-mer

        Returns:
            a tuple containing:
                - a read-only mapping from k-mer sequences to integer IDs
                - a read-only int64 array mapping each k-mer ID to its
                  reverse-complement ID
    """
    kmer_ids = {"".join(kmer): number
                for number, kmer in enumerate(product("ACGT",
                                                      repeat=kmer_size))}
    reverse_complement_ids = np.array(
        [kmer_ids[reverse_complement(kmer)] for kmer in kmer_ids],
        dtype=np.int64,
    )
    reverse_complement_ids.setflags(write=False)
    return MappingProxyType(kmer_ids), reverse_complement_ids


def _encode_kmers(promoter_seq: str,
                  kmer_size: int
                  ) -> tuple[np.ndarray, np.ndarray, np.ndarray]:
    """ Encodes overlapping k-mers, flags valid bases and encodes their
        reverse complements

        Each distinct k-mer of the given size is represented by an integer
        ID, with the same ID used across all promoters.

        Arguments:
            promoter_seq: the promoter DNA sequence to encode
            kmer_size: the number of bases in each k-mer

        Returns:
            a tuple of:
                - an int64 array of overlapping k-mer IDs, indexed by
                  their starting positions in the sequence
                - a bool array of the same length, True where the k-mer
                  contains only A, C, G or T
                - an int64 array of reverse-complement IDs in the same
                  order as the original k-mers

    """
    num_kmers = max(0, len(promoter_seq) - kmer_size + 1)
    encoded_kmers = np.zeros(num_kmers, dtype=np.int64)
    valid_kmers = np.ones(num_kmers, dtype=bool)
    if not num_kmers:
        return encoded_kmers, valid_kmers, np.zeros(0, dtype=np.int64)

    # get lookup - key: kmer_seq, value: kmer_id
    kmer_to_ids_dict, reverse_complement_ids = _get_kmer_to_ids(kmer_size)

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

        encoded_kmers[start] = current_kmer_id
    encoded_complements = reverse_complement_ids[encoded_kmers]
    return encoded_kmers, valid_kmers, encoded_complements


def count_repeat_motifs(promoter_seq: str,  # pylint: disable=too-many-locals
                        capped_max_count: int | None = config.MAX_MOTIFS_PER_PROMOTER,
                        *, bases_per_arm: int = config.BASES_PER_ARM,
                        max_spacer: int = config.MAX_SPACER
                        ) -> MotifCounts:
    """ Counts the inverted and direct repeat motifs in the given sequence

        A repeat is two arms separated by a spacer of 0 to max_spacer bases.
        In an inverted repeat the second arm is the reverse complement of the
        first, in a direct repeat the two arms are
        identical and get canonicalised. The spacer bases are not examined,
        so ACGNNCGT counts as an inverted repeat. Occurrences with an unknown
        base in either arm are skipped.

        Arguments:
            promoter_seq: the promoter sequence to scan, read 5' to 3'
            capped_max_count: the maximum count kept for each motif in this
                              promoter, or None to keep uncapped counts
            bases_per_arm: the number of bases in each repeat arm
            max_spacer: the maximum number of bases between repeat arms

        Returns:
            a MotifCounts instance with the sorted unique motif IDs and their
            counts, optionally capped at capped_max_count
    """

    # Skip promoters shorter than two adjacent repeat arms.
    if len(promoter_seq) < 2 * bases_per_arm:
        return MotifCounts.from_counts(np.zeros(0, dtype=np.int64),
                                       capped_max_count=capped_max_count)

    encoded_kmers, valid_kmers, encoded_complements = _encode_kmers(
        promoter_seq, bases_per_arm)

    # Each table stores counts at [arm_id, gap_size].
    inverted_counts = np.zeros((4 ** bases_per_arm, max_spacer + 1), dtype=np.int64)
    direct_counts = np.zeros_like(inverted_counts)

    # Check every possible arm pair for each allowed gap size.
    for gap_size in range(max_spacer + 1):
        second_arm_start_pos = bases_per_arm + gap_size
        first_arm_ids = encoded_kmers[:-second_arm_start_pos]

        # Keep only positions where both arms have known bases
        both_arms_valid = (
            valid_kmers[:-second_arm_start_pos]
            & valid_kmers[second_arm_start_pos:]
        )

        # Find inverted repeats
        is_inverted = first_arm_ids == encoded_complements[second_arm_start_pos:]
        inverted_arm_ids = first_arm_ids[both_arms_valid & is_inverted]

        # Find direct repeats (both arms are identical)
        is_direct = both_arms_valid & (
            first_arm_ids == encoded_kmers[second_arm_start_pos:])
        direct_arm_ids = np.minimum(
            first_arm_ids[is_direct],
            encoded_complements[:-second_arm_start_pos][is_direct],
        )
        # Count occurrences of each arm for this gap size.
        inverted_counts[:, gap_size] = np.bincount(inverted_arm_ids,
                                                   minlength=inverted_counts.shape[0])
        direct_counts[:, gap_size] = np.bincount(direct_arm_ids,
                                                 minlength=direct_counts.shape[0])

    # Flatten each arm's gap counts in order, with inverted repeats before
    # direct repeats. This preserves the motif IDs used by the scoring pipeline
    return MotifCounts.from_counts(
        np.concatenate((inverted_counts.ravel(), direct_counts.ravel())),
        capped_max_count=capped_max_count,
    )


def count_linear_motifs(promoter_seq: str,
                        capped_max_count: int | None = config.MAX_MOTIFS_PER_PROMOTER,
                        *, linear_length: int = config.LINEAR_LENGTH,
                        motif_id_start: int = config.linear_motif_id_start
                        ) -> MotifCounts:
    """ Counts the linear motifs of the given length in the sequence

        A motif and its reverse complement are counted as one motif,
        represented by the sequence with the lower k-mer ID.
        K-mers containing bases other than A, C, G or T are excluded from the counts.

        Arguments:
            promoter_seq: the promoter sequence to scan
            capped_max_count: the maximum count kept for each motif in this
                              promoter, or None to keep uncapped counts
            linear_length: the number of bases in each linear motif
            motif_id_start: the first ID reserved for linear motifs

        Returns:
            a MotifCounts instance with the sorted unique motif IDs and their
            counts, optionally capped at capped_max_count
    """

    if len(promoter_seq) < linear_length:
        return MotifCounts.from_counts(np.zeros(0, dtype=np.int64),
                                       capped_max_count=capped_max_count)
    encoded_kmers, valid_kmers, encoded_complements = _encode_kmers(
        promoter_seq, linear_length)
    # Use the smaller ID of each k-mer/reverse-complement pair so both
    # contribute to the same motif count.
    canonical_ids = np.minimum(encoded_kmers, encoded_complements)
    counts = np.bincount(canonical_ids[valid_kmers])
    motif_counts = MotifCounts.from_counts(counts, capped_max_count=capped_max_count)
    motif_counts.motif_ids += motif_id_start
    return motif_counts
