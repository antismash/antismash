from typing import List, Tuple, Dict, Any

from .data_structures import Prediction
from antismash.common.html_renderer import Markup
from antismash.detection.nrps_pks_domains import ModularDomain
from paras.run import run_paras as paras


class ParasResult(Prediction):
    """ Holds all the relevant results from PARAS for an adenylation domain """
    def __init__(self, paras_predictions: List[Tuple[float, str]]) -> None:
        super().__init__("PARAS")
        self.predictions = paras_predictions
        self.best_prediction = self.get_best_prediction()

    def get_classification(self) -> List[str]:
        if self.predictions:
            predictions = [self.predictions[0]]
            for prediction in self.predictions[1:]:
                if prediction[0] == predictions[0][0]:
                    predictions.append(prediction)
            return [prediction[1] for prediction in predictions]

        else:
            return []

    def get_best_prediction(self) -> str:
        if self.predictions:
            best_prediction = self.predictions[0][1]
        else:
            best_prediction = ''
        return best_prediction

    def as_html(self) -> Markup:
        raw_start = (
            "\n"
            "<dl><dt>PARAS prediction, score (0-1):</dt>\n"
            " <dd>\n"
            "  <dl>\n"
        )
        core = "\n".join("  <dd></dd><dt>%s: %.2f</dt>\n" % (name, score) for score, name in self.predictions)
        raw_end = (
            "\n"
            "  </dl>\n"
            " </dd>\n"
            "</dl>\n"
        )
        return Markup("%s%s%s" % (raw_start, core, raw_end))

    def to_json(self) -> Dict[str, Any]:
        return vars(self)

    @classmethod
    def from_json(cls, json: Dict[str, Any]) -> "ParasResult":
        """ Creates a Prediction from a JSON representation """
        return ParasResult(json["predictions"])


def run_paras(a_domains: List[ModularDomain]) -> Dict[str, "ParasResult"]:
    results: Dict[str, ParasResult] = {}
    for domain in a_domains:
        predictions = paras(domain.translation, threshold=0.5)
        results[domain.domain_id] = ParasResult(predictions)

    return results
