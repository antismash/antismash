# License: GNU Affero General Public License v3 or later
# A copy of GNU AGPL v3 should have been included in this software package in LICENSE.txt.

import logging
from typing import List, Set, Tuple

import Bio.Data.IUPACData
from Bio.SeqUtils.ProtParam import ProteinAnalysis
from Bio.SeqFeature import FeatureLocation

from .secmet import Feature, Record


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


def distance_to_pfam(record: Record, query: Feature, hmmer_profiles: List[str]) -> int:
    """ Checks how many nucleotides a gene is away from another gene with one
        of the given Pfams.

        Arguments:
            record: the secmet.Record containing the query and other genes to compare
            query: a secmet.Feature to use the location of
            hmmer_profiles: a list of profile names to search for in CDSFeature
                            domain lists

        Returns:
            the distance (in bases) to the closest gene with one of the profiles
            or -1 if no gene matching was found
    """
    max_range = 40000
    # Get all CDS features in record
    cds_features = record.get_cds_features()
    # Get all CDS features within <X nt distances
    close_cds_features = []
    distance = {}
    # TODO: change from O(n) to binary search for start position
    for cds in cds_features:
        search_range = FeatureLocation(query.location.start - max_range,
                                       query.location.end - max_range)
        if cds.is_contained_by(search_range):
            close_cds_features.append(cds)
            distance[cds.get_name()] = min([
                                abs(cds.location.start - query.location.end),
                                abs(cds.location.end - query.location.start),
                                abs(cds.location.start - query.location.start),
                                abs(cds.location.end - query.location.end)])
    # For nearby CDS features, check if they have hits to the pHMM
    closest_distance = -1
    for cds in close_cds_features:
        if cds.sec_met:
            for profile in hmmer_profiles:
                if profile in cds.sec_met.domains:
                    if closest_distance == -1 or distance[cds.get_name()] < closest_distance:
                        closest_distance = distance[cds.get_name()]
    return closest_distance
