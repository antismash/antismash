# License: GNU Affero General Public License v3 or later
# A copy of GNU AGPL v3 should have been included in this software package in LICENSE.txt.

""" Provides analysis of AT domain signatures """

from typing import Any, Dict, List, Tuple

from jinja2 import Markup

from antismash.common import path, subprocessing, utils, fasta
from antismash.modules.nrps_pks.data_structures import Prediction
from antismash.modules.nrps_pks.pks_names import get_long_form

_KS_REFERENCE_ALIGNMENT = path.get_full_path(__file__, "data", "649KS_sequences_031218.fasta")
_LEAF2CLADE_TBL = path.get_full_path(__file__, "data", "transPACT_leaf2clade.tsv")

class KSResult:
    """ A result for a specific KS domain """
    __slots__ = ["name", "clade", "specificity"]

    def __init__(self, name: str, clade: str, specificity: str) -> None:
        assert isinstance(name, str)
        assert isinstance(clade, str)
        assert isinstance(specificity, str)
        self.name = name
        self.clade = clade
        self.specificity = specificity

    def __repr__(self) -> str:
        return str(self)

    def __str__(self) -> str:
        return "KSResult(name=%s, clade=%s, specificity=%s)" % (self.name, self.clade, self.specificity)

    def to_json(self) -> Tuple[str, str, str]:
        """ Serialises the instance """
        return (self.name, self.clade, self.specificity)

    @staticmethod
    def from_json(json: Tuple[str, str, str]) -> "KSResult":
        """ Deserialise an KSResult instance """
        assert len(json) == 3
        return KSResult(*json)


class KSPrediction(Prediction):
    """ Holds the transPACT predictions for a domain"""
    def __init__(self, prediction: str) -> None:
        super().__init__("transPACT_KS")
        self.prediction = prediction
        
    def as_html(self) -> Markup:
        if not self.prediction:
            return Markup("No matches")
        line = "<dd>%s</dd>\n" % (self.prediction)
        html = ((
            "<dl>\n"
            " <dt>transPACT assigned specificiy:</dt>\n"
            "%s"
            "</dl>\n"
        ) % line)
        return Markup(html)

    def to_json(self) -> Dict[str, Any]:
        return {
            "method": "transPACT_KS",
            "prediction": self.prediction,
        }

    @staticmethod
    def from_json(json: Dict[str, Any]) -> "KSPrediction":
        assert json["method"] == "transPACT_KS"
        return KSPrediction(json["prediction"])


def run_transpact_ks_analysis(domains: Dict[str, str]) -> Dict[str, Prediction]:
    """ Analyses PKS signature of KS domains

        Arguments:
            domains: a dictionary mapping domain identifier (e.g. 'locus_KS2')
                     to domain sequence

        Returns:
            a dictionary mapping domain identifier to
                a list of KSResults
    """

    print("transPACT goes here!")
    pass
