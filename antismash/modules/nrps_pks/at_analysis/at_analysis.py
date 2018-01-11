# License: GNU Affero General Public License v3 or later
# A copy of GNU AGPL v3 should have been included in this software package in LICENSE.txt.

""" Provides analysis of AT domain signatures """

from typing import Dict, List, Tuple

from antismash.common import path, subprocessing, utils, fasta

_SIGNATURE_LENGTH = 24
_AT_DOMAINS_FILENAME = path.get_full_path(__file__, "data", "AT_domains_muscle.fasta")
_AT_POSITIONS_FILENAME = path.get_full_path(__file__, "data", "ATpositions.txt")
_SIGNATURES_FILENAME = path.get_full_path(__file__, "data", "pks_signatures.fasta")
_REF_SEQUENCE = "P0AAI9_AT1"


class ATSignatureResults(dict):
    """ Holds a list of ATResult for each gene name """
    def to_json(self) -> Dict[str, List[Tuple[str, str, float]]]:
        """ Serialises the instance """
        results = {}
        for key, value in self.items():
            results[key] = [result.to_json() for result in value]
        return results

    @staticmethod
    def from_json(json) -> "ATSignatureResults":
        """ Deserialises an ATSignatureResults instance """
        results = ATSignatureResults()
        for key, value in json.items():
            results[key] = [ATResult.from_json(val) for val in value]
        return results


class ATResult:
    """ A result for a specific AT domain """
    __slots__ = ["name", "signature", "score"]

    def __init__(self, name, signature, score):
        assert isinstance(name, str)
        assert isinstance(signature, str)
        assert isinstance(score, float)
        self.name = name
        self.signature = signature
        self.score = score

    def __repr__(self) -> str:
        return str(self)

    def __str__(self) -> str:
        return "ATResult(query=%s, signature=%s, score=%.1f)" % (self.name, self.signature, self.score)

    def to_json(self) -> Tuple[str, str, float]:
        """ Serialises the instance """
        return (self.name, self.signature, self.score)

    @staticmethod
    def from_json(json) -> "ATResult":
        """ Deserialise an ATResult instance """
        return ATResult(*json)


def get_at_positions(startpos=7):
    """ Reads a reference list of positions used in signature extraction
        from file.
    """
    with open(_AT_POSITIONS_FILENAME, "r") as handle:
        text = handle.read().strip().replace(' ', '_')
    positions = [int(pos) - startpos for pos in text.split("\t")]
    return positions


def score_signatures(query_signatures: Dict[str, str],
                     reference_signatures: Dict[str, str]
                     ) -> Dict[str, List[ATResult]]:
    """ Scores PKS signature by comparing against database of signatures.

        The score is calculated as a percentage of pairwise matches for the
        length of the sequence. Only scores greater than 50% are considered.

        Arguments:
            query_signatures: a dictionary mapping query name to signature
            reference_signatures: a dictionary mapping reference name to signature

        Returns:
            a wrapped dictionary mapping query name to at most 10 ATResult instances,
            sorted in order of decreasing score
    """
    results = ATSignatureResults()
    for key, query_sig_seq in sorted(query_signatures.items()):
        scores = []
        for sig_name, sig_seq in reference_signatures.items():
            score = 0.
            for query_amino, sig_amino in zip(query_sig_seq, sig_seq):
                if query_amino == sig_amino:
                    score += 1.
            score = 100 * score / _SIGNATURE_LENGTH
            # ignore scores <= 50%
            if score > 50.:
                scores.append(ATResult(sig_name, reference_signatures[sig_name], score))
        # limit to 10 best hits, scores descending, names ascending for ties
        results[key] = sorted(scores, key=lambda x: x.score, reverse=True)[:10]
    return results


def run_at_domain_analysis(domains: Dict[str, str]) -> ATSignatureResults:
    """ Analyses PKS signature of AT domains

        Arguments:
            domains: a dictionary mapping domain identifier (e.g. 'locus_AT2')
                     to domain sequence

        Returns:
            a dictionary mapping domain identifier to
                a list of ATResults ordered by decreasing score
    """
    # construct the query signatures
    query_signatures = {}
    at_positions = get_at_positions(startpos=7)
    for name, seq in sorted(domains.items()):
        alignments = subprocessing.run_muscle_single(name, seq, _AT_DOMAINS_FILENAME)
        query_signatures[name] = utils.extract_by_reference_positions(alignments[name],
                                         alignments[_REF_SEQUENCE], at_positions)
    # load reference PKS signatures and score queries against them
    return score_signatures(query_signatures, fasta.read_fasta(_SIGNATURES_FILENAME))
