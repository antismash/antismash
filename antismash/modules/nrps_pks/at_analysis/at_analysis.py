# License: GNU Affero General Public License v3 or later
# A copy of GNU AGPL v3 should have been included in this software package in LICENSE.txt.

""" Provides analysis of AT domain signatures """

from typing import Any, Dict, List, Tuple

from jinja2 import Markup

from antismash.common import path, subprocessing, utils, fasta
from antismash.modules.nrps_pks.data_structures import Prediction
from antismash.modules.nrps_pks.pks_names import get_long_form

_SIGNATURE_LENGTH = 24
_AT_DOMAINS_FILENAME = path.get_full_path(__file__, "data", "AT_domains_muscle.fasta")
_AT_POSITIONS_FILENAME = path.get_full_path(__file__, "data", "ATpositions.txt")
_SIGNATURES_FILENAME = path.get_full_path(__file__, "data", "pks_signatures.fasta")
_REF_SEQUENCE = "P0AAI9_AT1"


class ATResult:
    """ A result for a specific AT domain """
    __slots__ = ["name", "signature", "score"]

    def __init__(self, name: str, signature: str, score: float) -> None:
        assert isinstance(name, str)
        assert isinstance(signature, str)
        assert isinstance(score, float)
        self.name = name
        self.signature = signature
        self.score = score

    def __repr__(self) -> str:
        return str(self)

    def __str__(self) -> str:
        return "ATResult(name=%s, signature=%s, score=%.1f)" % (self.name, self.signature, self.score)

    def to_json(self) -> Tuple[str, str, float]:
        """ Serialises the instance """
        return (self.name, self.signature, self.score)

    @staticmethod
    def from_json(json: Tuple[str, str, float]) -> "ATResult":
        """ Deserialise an ATResult instance """
        assert len(json) == 3
        return ATResult(*json)


class ATPrediction(Prediction):
    """ Holds the signature-based predictions for a domain"""
    def __init__(self, predictions: Dict[str, ATResult]) -> None:
        super().__init__("ATSignature")
        self.predictions = sorted(predictions.items(), key=lambda x: (-x[1].score, x[1].name))

    def get_classification(self) -> List[str]:
        results = []  # type: List[str]
        if not self.predictions:
            return results
        best_score = self.predictions[0][1].score
        for monomer, pred in self.predictions:
            if pred.score < best_score:
                break
            results.append(get_long_form(monomer))
        return results

    def as_html(self) -> Markup:
        if not self.predictions:
            return Markup("No matches above 50%")
        lines = []
        for monomer, pred in self.predictions[:3]:
            lines.append("<dd>%s: %.1f%%</dd>\n" % (get_long_form(monomer), pred.score))
        html = ((
            "<dl>\n"
            " <dt>Top 3 matches:</dt>\n"
            "%s"
            "</dl>\n"
        ) % "".join(lines))
        return Markup(html)

    def to_json(self) -> Dict[str, Any]:
        return {
            "method": "ATSignature",
            "predictions": {monomer: pred.to_json() for monomer, pred in self.predictions},
        }

    @staticmethod
    def from_json(json: Dict[str, Any]) -> "ATPrediction":
        assert json["method"] == "ATSignature"
        return ATPrediction({monomer: ATResult.from_json(pred) for monomer, pred in json["predictions"].items()})


def get_at_positions(startpos: int = 7) -> List[int]:
    """ Reads a reference list of positions used in signature extraction
        from file.
    """
    with open(_AT_POSITIONS_FILENAME, "r") as handle:
        text = handle.read().strip().replace(' ', '_')
    positions = [int(pos) - startpos for pos in text.split("\t")]
    return positions


def score_signatures(query_signatures: Dict[str, str],
                     reference_signatures: Dict[str, str]) -> Dict[str, Prediction]:
    """ Scores PKS signature by comparing against database of signatures.

        The score is calculated as a percentage of pairwise matches for the
        length of the sequence. Only scores greater than 50% are considered.

        Arguments:
            query_signatures: a dictionary mapping query name to signature
            reference_signatures: a dictionary mapping reference name to signature

        Returns:
            a dictionary mapping each query identifier to an ATPrediction
    """
    results = {}  # type: Dict[str, Prediction]
    for key, query_sig_seq in sorted(query_signatures.items()):
        # keep a single best prediction for each monomer type
        scores = {}  # type: Dict[str, ATResult]
        for sig_name, sig_seq in reference_signatures.items():
            monomer = sig_name.rsplit('_', 1)[-1]
            score = 0.
            for query_amino, sig_amino in zip(query_sig_seq, sig_seq):
                if query_amino == sig_amino:
                    score += 1.
            score = 100 * score / _SIGNATURE_LENGTH
            # ignore scores <= 50%
            if score > (scores[monomer].score if monomer in scores else 50.):
                scores[monomer] = ATResult(sig_name, reference_signatures[sig_name], score)
        results[key] = ATPrediction(scores)
    return results


def run_at_domain_analysis(domains: Dict[str, str]) -> Dict[str, Prediction]:
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
