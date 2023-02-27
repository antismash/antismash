from typing import List, Dict, Any

from .data_structures import Prediction
from antismash.common.html_renderer import Markup
from antismash.detection.nrps_pks_domains import ModularDomain
from paras.run import run_paras_bulk as paras
from paras.features import get_metadata, ParasSubstrate
from dataclasses import dataclass

GOOD_PARAS_THRESHOLD = 0.6
MINIMUM_PARAS_THRESHOLD = 0.2


@dataclass
class ParasPrediction:
    full_name: str
    paras_abbreviation: str
    norine_abbreviation: str
    probability: float

    @classmethod
    def from_paras_substrate(cls, paras_substrate: ParasSubstrate, probability: float) -> "ParasPrediction":
        return cls(paras_substrate.name, paras_substrate.abbreviation, paras_substrate.norine, probability)

    def to_json(self):
        return dict(vars(self))


class ParasResult(Prediction):
    """ Holds all the relevant results from PARAS for an adenylation domain """
    def __init__(self, paras_predictions: List[ParasPrediction]) -> None:
        super().__init__("PARAS")
        self.predictions = paras_predictions
        self.predictions.sort(key=lambda x: x.probability, reverse=True)

    def get_classification(self, as_norine: bool = False) -> List[str]:
        if self.predictions and self.predictions[0].probability >= GOOD_PARAS_THRESHOLD:
            predictions = [self.predictions[0]]
            for prediction in self.predictions[1:]:
                if prediction.probability == predictions[0].probability:
                    predictions.append(prediction)
            print([prediction.norine_abbreviation for prediction in predictions])
            return [prediction.norine_abbreviation for prediction in predictions]

        else:
            return []

    def as_html(self) -> Markup:
        if not self.predictions:
            return Markup("No hits above threshold.")

        raw_start = (
            "\n"
            "<dl><dt>Prediction probabilities:</dt>\n"
            " <dd><dl>\n"
        )
        core = "\n".join("  <dd></dd><dt>%s: %.2f</dt>\n" % (prediction.full_name, prediction.probability) for prediction in self.predictions)
        raw_end = (
            "  </dl></dd></dl>"
        )
        return Markup("%s%s%s" % (raw_start, core, raw_end))

    def to_json(self) -> Dict[str, Any]:
        return {"predictions": [dict(vars(prediction)) for prediction in self.predictions]}

    @classmethod
    def from_json(cls, json: Dict[str, Any]) -> "ParasResult":
        """ Creates a Prediction from a JSON representation """
        return cls(json["predictions"])


def run_paras(a_domains: List[ModularDomain]) -> Dict[str, Prediction]:
    results: Dict[str, Prediction] = {}
    sequences = []
    for domain in a_domains:
        sequences.append(domain.translation)
    predictions = paras(sequences, threshold=MINIMUM_PARAS_THRESHOLD)
    for i, domain in enumerate(a_domains):
        prediction = predictions[i]
        if domain.domain_id is not None:
            paras_predictions = []
            for probability, substrate in prediction:
                paras_substrate = get_metadata(substrate)
                paras_prediction = ParasPrediction.from_paras_substrate(paras_substrate, probability)
                paras_predictions.append(paras_prediction)

            paras_result = ParasResult(paras_predictions)
            results[domain.domain_id] = paras_result

    return results