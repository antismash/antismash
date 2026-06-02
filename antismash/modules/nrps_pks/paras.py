# License: GNU Affero General Public License v3 or later
# A copy of GNU AGPL v3 should have been included in this software package in LICENSE.txt.

""" Predicts A domain specificity using PARAS and PARASECT models """

import os
from typing import Any
from enum import Enum

from joblib import load

from parasect.core.retrain_models import retrain_model, model_needs_retraining, update_metadata_file
from parasect.core.models import ModelType
from parasect.core.parasect_result import MinimalResult
from parasect.core.antismash import run_parasect_minimal, run_paras_minimal

from antismash.config import ConfigType, get_config
from antismash.common.path import find_latest_database_version
from antismash.common.html_renderer import Markup
from antismash.detection.nrps_pks_domains import ModularDomain

from .data_structures import Prediction
from .signatures import get_a_dom_signatures
from .name_mappings import get_substrate_by_name, SubstrateName

PARASECT_CONFIDENCE_THRESHOLD = 0.7
PARAS_CONFIDENCE_THRESHOLD = 0.5


def _get_model_dir(config: ConfigType) -> str:
    """ A helper to construct the absolute path to the NRPS PARAS model base dir in the
        data directory.

        Arguments:
            config: the antiSMASH config

        Returns:
            the absolute path of the NRPS PARAS model base dir
    """
    root = os.path.join(config.database_dir, "nrps_pks", "paras")
    version = find_latest_database_version(root)
    return os.path.join(root, version)


def prepare_data(logging_only: bool = False) -> list[str]:
    """ Ensures packaged data is fully prepared

        Arguments:
            options: antiSMASH configuration options
            logging_only: whether to return error messages instead of raising exceptions

        Returns:
            a list of error messages (only if logging_only is True)
    """
    options = get_config()
    failure_messages: list[str] = []

    if "mounted_at_runtime" in options.database_dir:
        return failure_messages

    try:
        model_dir = _get_model_dir(options)
        metadata_path = os.path.join(model_dir, "metadata.txt")

        if not os.path.exists(metadata_path):
            if not logging_only:
                raise FileNotFoundError
            failure_messages.append(f"Failed to locate {metadata_path}")
        else:
            retrain_paras_models_if_needed(metadata_path, model_dir)

    except (FileNotFoundError, ValueError):
        if not logging_only:
            raise
        failure_messages.append("Failed to locate PARAS model dir.")

    return failure_messages


def check_prereqs(options: ConfigType) -> list[str]:
    """ Check if the metadata file is present. """
    failure_messages: list[str] = []
    if "mounted_at_runtime" in options.database_dir:  # can't prepare this one
        return failure_messages
    try:
        model_dir = _get_model_dir(options)
    except ValueError:
        model_dir = ""
    if not os.path.exists(model_dir):
        failure_messages.append(f"Failed to locate {model_dir}")
    return failure_messages


def retrain_paras_models_if_needed(metadata_path: str, model_dir: str) -> None:
    """ Retrains PARAS models if sklearn versions do not match
    Args:
        metadata_path: path to PARAS model metadata file
        model_dir: path to PARAS model dir
    """
    for model_type in ModelType.ANTISMASH_MODELS:
        if model_needs_retraining(metadata_path, model_type):
            model = retrain_model(model_type)
            model.save(model_dir)
            update_metadata_file(model_type, metadata_path)


class ParasModel(Enum):
    """ Enum for PARAS model types"""
    PARAS = 1
    PARASECT = 2


class GeneralParasResult(Prediction):
    """ Holds all the relevant results from PARASECT for a domain """

    def __init__(self, predictor: ParasModel,
                 aa34: str,
                 substrate_predictions: list[str],
                 prediction_confidences: list[float],
                 substrates_to_report: int = 3) -> None:
        super().__init__(predictor.name.lower())
        assert len(substrate_predictions) == len(prediction_confidences)
        assert len(substrate_predictions) >= substrates_to_report

        self.aa34 = aa34
        self.predicted_substrate = substrate_predictions[0]
        self.confidence = prediction_confidences[0]
        self.top_substrates = substrate_predictions[:substrates_to_report]
        self.top_confidences = prediction_confidences[:substrates_to_report]

        self.uncertain = False

        if predictor == ParasModel.PARAS and self.confidence < PARAS_CONFIDENCE_THRESHOLD:
            self.uncertain = True
        if predictor == ParasModel.PARASECT and self.confidence < PARASECT_CONFIDENCE_THRESHOLD:
            self.uncertain = True

    def _get_classification(self) -> list[SubstrateName]:
        classification: list[SubstrateName] = []

        for i, substrate in enumerate(self.top_substrates):
            confidence = self.top_confidences[i]
            if self.confidence - confidence < 0.2:
                classification.append(get_substrate_by_name(substrate))


        return classification

    def get_classification(self, as_norine: bool = False) -> list[str]:
        """ Get the classification """
        substrates = self._get_classification()

        def mapper(substrate: SubstrateName) -> str:
            if as_norine:
                return substrate.norine
            return substrate.norine if substrate.norine != "X" else substrate.short

        return list(map(mapper, substrates))

    def as_html(self) -> Markup:
        note = ""
        if self.uncertain:
            note = "<strong>NOTE: low confidence prediction</strong><br>\n"

        raw = ("\n"
               f"<dl><dt>{self.method.upper()} details:</dt>\n"
               " <dd>"
               f"  {note}"
               "  <dl>"
               "   <dt>Predicted substrate (confidence):</dt>\n")

        for i, substrate in enumerate(self.top_substrates):
            confidence = self.top_confidences[i]
            raw += f"   <dd>{substrate} ({confidence * 100:.1f}%)</dd>\n"

        raw += ("  </dl>\n"
                " </dd>\n"
                "</dl>\n")


        return Markup(raw)

    def __str__(self) -> str:
        return f"{self.method.capitalize()}Result: " + str(vars(self))

    def to_json(self) -> dict[str, Any]:
        return {
            "aa34": self.aa34,
            "predicted_substrate": self.predicted_substrate,
            "confidence": self.confidence,
            "top_substrates": self.top_substrates,
            "top_confidences": self.top_confidences
        }

    @classmethod
    def from_json(cls, json: dict[str, Any]) -> "Prediction":
        """ Creates a Prediction from a JSON representation """
        raise NotImplementedError(f"Prediction subclass {cls} did not implement from_json()")

class ParasResult(GeneralParasResult):
    """ Holds all the relevant results from PARAS for a domain """

    def __init__(self,
                 aa34: str,
                 substrate_predictions: list[str],
                 prediction_confidences: list[float],
                 substrates_to_report: int = 1) -> None:
        super().__init__(ParasModel.PARAS, aa34, substrate_predictions,
                         prediction_confidences, substrates_to_report)

    @classmethod
    def from_paras_result(cls, minimal_result: MinimalResult,
                          substrates_to_report: int = 1) -> "ParasResult":
        """ Generates a PARAS Result from a minimal result"""

        return cls(minimal_result.signature, minimal_result.prediction_labels,
                   minimal_result.predictions,
                   substrates_to_report=substrates_to_report)

    @classmethod
    def from_json(cls, json: dict[str, Any]) -> "ParasResult":

        return ParasResult(json["aa34"],
                           json["top_substrates"],
                           json["top_confidences"],
                           len(json["top_substrates"]))


class ParasectResult(GeneralParasResult):
    """ Holds all the relevant results from PARASECT for a domain """
    def __init__(self, aa34: str,
                 substrate_predictions: list[str],
                 prediction_confidences: list[float],
                 substrates_to_report: int = 3):
        super().__init__(ParasModel.PARASECT, aa34, substrate_predictions,
                         prediction_confidences, substrates_to_report)

    @classmethod
    def from_parasect_result(cls, minimal_result: MinimalResult,
                          substrates_to_report: int = 3) -> "ParasectResult":
        """ Generates a PARASECT Result from a minimal result"""

        return cls(minimal_result.signature, minimal_result.prediction_labels,
                   minimal_result.predictions,
                   substrates_to_report=substrates_to_report)

    @classmethod
    def from_json(cls, json: dict[str, Any]) -> "ParasectResult":
        return ParasectResult(json["aa34"],
                              json["top_substrates"],
                              json["top_confidences"],
                              len(json["top_substrates"]))

def run_paras(a_domains: list[ModularDomain], options: ConfigType) -> dict[str, Prediction]:
    """ Runs nrpys over the provided A domains.

        Arguments:
            a_domains: a list of ModularDomains, one for each A domain
            options: antismash options

        Returns:
            a dictionary mapping each domain name to a PredictorSVMResult
    """

    model_path = os.path.join(_get_model_dir(options), "all_substrates_model.paras.gz")
    model = load(model_path)

    names: list[str] = []
    signatures: list[str] = []
    for sig, domain in [(get_a_dom_signatures(a_domain)[1], a_domain) for a_domain in a_domains]:
        if not sig:
            continue

        names.append(domain.get_name())
        signatures.append(sig)

    paras_results = run_paras_minimal(model, signatures, names)

    results: dict[str, Prediction] = {}

    for paras_result in paras_results:
        result = ParasResult.from_paras_result(paras_result)
        results[paras_result.domain_name] = result

    return results


def run_parasect(a_domains: list[ModularDomain], options: ConfigType) -> dict[str, Prediction]:
    """ Runs parasect over the provided A domains.

        Arguments:
            a_domains: a list of ModularDomains, one for each A domain
            options: antismash options

        Returns:
            a dictionary mapping each domain name to a PredictorSVMResult
    """

    model_path = os.path.join(_get_model_dir(options), "bacterial_model.parasect.gz")
    model = load(model_path)

    names: list[str] = []
    signatures: list[str] = []

    for sig, domain in [(get_a_dom_signatures(a_domain)[1], a_domain) for a_domain in a_domains]:
        if not sig:
            continue

        names.append(domain.get_name())
        signatures.append(sig)

    parasect_results = run_parasect_minimal(model, signatures, names, bacterial_only=True)

    results: dict[str, Prediction] = {}

    for parasect_result in parasect_results:
        result = ParasectResult.from_parasect_result(parasect_result)
        results[parasect_result.domain_name] = result

    return results
