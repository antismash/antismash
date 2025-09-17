# License: GNU Affero General Public License v3 or later
# A copy of GNU AGPL v3 should have been included in this software package in LICENSE.txt.

""" A collection of functions and classes that are used by multiple modules
    within antismash but aren't large enough to justify their own module within
    antismash.common.
"""

import dataclasses
from typing import Dict, Iterable, List, Optional

import Bio.Data.IUPACData
from Bio.SearchIO import QueryResult
from Bio.SeqUtils.ProtParam import ProteinAnalysis

from .fasta import read_fasta
from .secmet import Feature, Record


@dataclasses.dataclass(kw_only=True, slots=True, eq=True)
class Hit:
    """ A container for a pairwise match of some query to a reference within a database.
        E.g. a pHMM hit or a BLAST/diamond hit.
    """
    query_id: str
    reference_id: str
    identity: Optional[float] = None
    evalue: Optional[float] = None
    bitscore: Optional[float] = None


class RobustProteinAnalysis(ProteinAnalysis):
    """ A simple subclass of ProteinAnalysis that can deal with
        a protein sequence containing invalid characters.

        If ignoring invalid characters, the molecular weight is increased by
        the average weight of an amino-acid (i.e. 110) for each invalid case.
    """
    PROTEIN_LETTERS = set(Bio.Data.IUPACData.protein_letters)

    def __init__(self, prot_sequence: str, monoisotopic: bool = False,
                 ignore_invalid: bool = True) -> None:
        if not isinstance(ignore_invalid, bool):
            raise TypeError("ignore_invalid must be a boolean")
        self._ignore_invalid = ignore_invalid

        prot_sequence = prot_sequence.upper()

        self.original_sequence = prot_sequence
        # remove all invalids
        prot_sequence = "".join(filter(lambda x: x in RobustProteinAnalysis.PROTEIN_LETTERS,
                                       self.original_sequence))
        super().__init__(prot_sequence, monoisotopic)

    def molecular_weight(self, *, weighting_multiplier: float = 110) -> float:
        weight = super().molecular_weight()
        if not self._ignore_invalid:
            aa_difference = len(self.original_sequence) - len(self.sequence)
            weight += weighting_multiplier * aa_difference
        return weight


def extract_by_reference_positions(query: str, reference: str, ref_positions: List[int]) -> Optional[str]:
    """ Extracts the given positions from a query alignment. The positions are
        adjusted to account for any gaps in the reference sequence.

        Arguments:
            query: the aligned query
            reference: the aligned reference
            ref_positions: the positions of interest in the unaligned reference

        Returns:
            a string containing the sequence elements at the adjusted reference
            positions or None if the reference is too short for some reason
    """
    # adjust position of interest to account for gaps in the ref sequence alignment
    positions = []
    position_skipping_gaps = 0
    for i, amino in enumerate(reference):
        if amino in "-.":
            continue
        if position_skipping_gaps in ref_positions:
            positions.append(i)
        position_skipping_gaps += 1
    if len(positions) != len(ref_positions):
        return None
    # extract positions from query sequence
    return "".join([query[i] for i in positions])


def extract_from_alignment(hit: QueryResult, positions: Iterable[int]) -> Optional[str]:
    """ Extracts bases from the given alignment of a sequence against a profile,
        at the given positions. The positions are adjusted to account for any
        gaps in the reference sequence.

        If not all positions can be satisfied, no result will be returned.

        Arguments:
            hit: the hit with the alignment
            positions: the reference positions of the bases to extract

        Returns:
            a string containing the extracted bases
            or None if not all positions could be extracted
    """
    profile = hit.aln[1].seq
    query = hit.aln[0].seq
    offset = hit.hit_start
    return extract_by_reference_positions(query, profile, [p - offset for p in positions if p >= offset])


def distance_to_pfam(record: Record, query: Feature, hmmer_profiles: list[str],
                     *, max_range: int = 40_000,
                     ) -> int:
    """ Checks how many nucleotides a gene is away from another gene with one
        of the given Pfams.

        Arguments:
            record: the secmet.Record containing the query and other genes to compare
            query: a secmet.Feature to use the location of
            hmmer_profiles: a list of profile names to search for in CDSFeature
                            domain lists
            max_range: the maximum search range in nucleotides

        Returns:
            the distance (in bases) to the closest gene with one of the profiles
            or -1 if no gene matching was found
    """
    if record.is_circular():
        max_range = min(max_range, len(record) // 2)
    search_range = record.extend_location(query.location, max_range)
    close_cds_features = record.get_cds_features_within_location(search_range, with_overlapping=True)

    # For nearby CDS features, check if they have hits to the pHMM
    profiles = set(hmmer_profiles)
    closest_distance = -1
    for cds in close_cds_features:
        if cds.sec_met is None:
            continue
        if set(cds.sec_met.domain_ids).intersection(profiles):
            distance = record.get_distance_between_features(query, cds)
            if closest_distance == -1 or distance < closest_distance:
                closest_distance = distance
    return closest_distance


def get_hmm_lengths(hmm_file: str) -> Dict[str, int]:
    """ Finds the lengths of all HMM profiles in the provided HMM file.

        Arguments:
            hmm_file: the HMM file to read

        Returns:
            a dictionary mapping each NAME field in the file to it's LENG field
    """
    lengths = {}
    with open(hmm_file, "r", encoding="utf-8") as handle:
        contents = handle.read()
    contents = contents.replace("\r", "")
    hmms = contents.split("//")[:-1]
    for hmm in hmms:
        namepart = hmm.split("NAME  ")[1]
        name = namepart.split("\n")[0]
        lengthpart = hmm.split("LENG  ")[1]
        length = lengthpart.split("\n")[0]
        lengths[name] = int(length)
    return lengths


def get_fasta_lengths(fasta_file: str) -> Dict[str, int]:
    """ Finds the lengths of all sequences in the provided FASTA file.

        Arguments:
            hmm_file: the FASTA file to read

        Returns:
            a dictionary mapping sequence identifier to sequence length
    """
    fasta = read_fasta(fasta_file)
    return {key: len(val) for key, val in fasta.items()}
