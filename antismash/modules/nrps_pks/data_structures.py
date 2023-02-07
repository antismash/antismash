# License: GNU Affero General Public License v3 or later
# A copy of GNU AGPL v3 should have been included in this software package in LICENSE.txt.

""" Contains data structures for the nrps_pks module """

from typing import Any, Dict, List, Optional

from antismash.common.html_renderer import Markup


class Prediction:
    """ The base of all predictions for a domain. Must be provided with the method
        creating the prediction as a minimum. Can contain whatever extra information
        is relevant, but must be able to convert to and from JSON, along with providing
        a standardised DomainPrediction and HTML conversion.
    """
    def __init__(self, method: str) -> None:
        self.method = method

    def get_classification(self, as_norine: bool = False) -> List[str]:
        """ Returns a list of equally likely predictions. If no prediction could
            be made, an empty list is returned.
            The optional as_norine field forces the use of valid Norine names.
        """
        raise NotImplementedError(f"Prediction subclass {type(self)} did not implement get_classification()")

    def as_html(self) -> Markup:
        """ Returns a jinja2.Markup object containing HTML to use when representing
            this prediction in the sidepanel output.
        """
        raise NotImplementedError(f"Prediction subclass {type(self)} did not implement as_html()")

    def to_json(self) -> Dict[str, Any]:
        """ Creates a JSON representation of this prediction """
        raise NotImplementedError(f"Prediction subclass {type(self)} did not implement to_json()")

    @classmethod
    def from_json(cls, json: Dict[str, Any]) -> "Prediction":
        """ Creates a Prediction from a JSON representation """
        raise NotImplementedError(f"Prediction subclass {cls} did not implement from_json()")


class SimplePrediction(Prediction):
    """ The simplest case of a Prediction, containing a simple method"""
    def __init__(self, method: str, prediction: str,
                 norine_prediction: Optional[str] = None) -> None:
        super().__init__(method)
        self.prediction = prediction
        self.norine_prediction = norine_prediction if norine_prediction else prediction

    def get_classification(self, as_norine: bool = False) -> List[str]:
        if as_norine:
            return [self.norine_prediction]
        return [self.prediction]

    def as_html(self) -> Markup:
        return Markup(f"{self.method}: {self.prediction}")

    def to_json(self) -> Dict[str, Any]:
        return {"class": "SimplePrediction",
                "method": self.method,
                "prediction": self.prediction}

    @classmethod
    def from_json(cls, json: Dict[str, Any]) -> "SimplePrediction":
        assert json.get("class") == "SimplePrediction"
        return cls(json["method"], json["prediction"])
