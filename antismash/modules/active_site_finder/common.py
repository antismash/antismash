# License: GNU Affero General Public License v3 or later
# A copy of GNU AGPL v3 should have been included in this software package in LICENSE.txt.

""" Common functions and classes for use in active site finders various analyses """

import os
from typing import List, Optional, Tuple

from antismash.common import fasta, path, secmet, subprocessing, utils


class Alignment:  # pylint: disable=too-few-public-methods
    """ An analog for a hmmer HSP. """
    def __init__(self, domain: secmet.features.Domain, query: str,
                 profile: str, hit_start: int, hit_end: int) -> None:
        assert isinstance(domain, secmet.features.Domain)
        self.domain = domain
        self.query = str(query)
        self.profile = str(profile)
        self.hit_start = int(hit_start)
        self.hit_end = int(hit_end)
        if self.hit_start < 0:
            raise ValueError("Hit start position is negative")
        if self.hit_end <= self.hit_start:
            raise ValueError("Hit end position is not greater than the start position")

    def get_signature(self, positions: List[int]) -> str:
        """ Returns a signature extracted from the query sequence using the
            given positions within the reference profile.
            Positions should be 1-indexed.
        """
        return get_signature(self.query, self.profile, [position - self.hit_start for position in positions])


class ActiveSiteAnalysis:
    """ A generic analysis framework. """
    def __init__(self, target_domain: str, candidates: Tuple[secmet.features.Domain, ...],
                 database: str, positions: List[int],
                 expected_values: List[str], emissions: Optional[List[float]] = None) -> None:
        self.target_domain = str(target_domain)
        self.database = path.get_full_path(__file__, 'data', database)
        if not os.path.exists(self.database):
            raise ValueError(f"No database file located for: {self.database!r}")
        self.positions = list(map(int, positions))
        self.expected_values = list(map(str, expected_values))
        if len(expected_values) != len(positions):
            raise ValueError("Number of expected values must match number of positions")
        self.emissions = None
        if emissions:
            self.emissions = list(map(int, emissions))
            if len(self.emissions) != len(positions):
                raise ValueError("Number of emissions must match number of positions")

        self.domains_of_interest: List[secmet.features.Domain] = []
        for candidate in candidates:
            if not isinstance(candidate, secmet.features.Domain):
                raise TypeError(f"Candidates must be Domains, not {type(candidate)}")
            if candidate.domain == self.target_domain:
                self.domains_of_interest.append(candidate)

    def get_alignments(self) -> List[Alignment]:
        """ Builds an Alignment for each hit in the results of running the
            provided command on the provided data.
        """
        if not self.domains_of_interest:
            return []

        # for safety of the tools, rename long domain names to a simple numeric index
        data = fasta.get_fasta_from_features(self.domains_of_interest, numeric_names=True)
        assert data, "empty fasta created"

        extra_args = ["-T", "0",  # min score
                      "-E", "0.1"]  # max evalue
        results = subprocessing.run_hmmpfam2(self.database, data, extra_args=extra_args)

        alignments = []
        for result in results:
            if not result.hsps:
                continue
            assert result.id == result.hsps[0].aln[0].id
            # fetch back the real domain from the numeric index used in the fasta
            domain = self.domains_of_interest[int(result.id)]
            alignments.append(Alignment(domain, result.hsps[0].aln[0].seq, result.hsps[0].aln[1].seq,
                                        result.hsps[0].hit_start, result.hsps[0].hit_end))
        return alignments

    def scaffold_matches(self, alignment: Alignment) -> bool:
        """ Returns True if the query sequence in the given alignment matches the expected signature """
        return get_signature(alignment.query, alignment.profile, self.positions) == "".join(self.expected_values)


def get_signature(query: str, hmm: str, positions: List[int]) -> str:
    """ Retrieves a signature from an aligned pair based on 1-indexed positions given.

        Arguments:
            query: the sequence of the query that the signature will be extracted from
            hmm: the sequence of the hit, used to adjust positions to account for introduced gaps
            positions: a list of 1-indexed positions to use for the signature
                       positions are relative to hit start
            expected: if provided, a signature extracted from the reference sequence must
                      match this

        Returns:
            None if the provided positions would be out of bounds of the hit,
            otherwise a string of the same length as positions and expected
    """
    ungapped = str(hmm).replace('.', '')
    if max(positions) > len(ungapped):
        # the hit was too small and a correct signature can't be generated
        return ""
    signature = utils.extract_by_reference_positions(query, hmm, [pos - 1 for pos in positions])
    assert signature
    return signature
