from typing import List, Tuple, Dict, Any

from .data_structures import Prediction
from antismash.common.html_renderer import Markup
from antismash.detection.nrps_pks_domains import ModularDomain
from paras.run import run_paras_bulk


class ParasResult(Prediction):
    """ Holds all the relevant results from PARAS for an adenylation domain """
    def __init__(self, paras_predictions: List[Tuple[float, str]]) -> None:
        super().__init__("PARAS")
        self.predictions = paras_predictions

    def get_classification(self) -> List[str]:
        if not self.predictions:
            return []

        predictions = [self.predictions[0]]
        for prediction in self.predictions[1:]:
            if prediction[0] == predictions[0][0]:
                predictions.append(prediction)
        return [prediction[1] for prediction in predictions]

    def as_html(self) -> Markup:
        if not self.predictions:
            return Markup("No hits above threshold.")

        raw_start = (
            "<dl><dt>PARAS prediction, score (0-1):</dt>"
            " <dd>"
            "  <dl>"
        )
        core = "\n".join(f"<dd></dd><dt>{name}: {score:.2f}</dt>" for score, name in self.predictions)
        raw_end = (
            "  </dl>"
            " </dd>"
            "</dl>"
        )
        return Markup(f"{raw_start}{core}{raw_end}")

    def to_json(self) -> Dict[str, Any]:
        return vars(self)

    @classmethod
    def from_json(cls, json: Dict[str, Any]) -> Prediction:
        """ Creates a Prediction from a JSON representation """
        return ParasResult(json["predictions"])


def run_paras(a_domains: List[ModularDomain]) -> Dict[str, Prediction]:
    results: Dict[str, Prediction] = {}
    sequences = []
    for domain in a_domains:
        sequences.append(domain.translation)
    predictions = run_paras_bulk(sequences, threshold=0.2)
    for i, domain in enumerate(a_domains):
        assert domain.domain_id
        prediction = predictions[i]
        results[domain.domain_id] = ParasResult(prediction)

    return results
