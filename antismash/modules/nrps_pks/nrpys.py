# License: GNU Affero General Public License v3 or later
# A copy of GNU AGPL v3 should have been included in this software package in LICENSE.txt.

""" Provides a collection of functions and classes to run the external Rust program
    nrps_rs and interpret the results
"""

from dataclasses import dataclass
import itertools
import os
from typing import Any, Optional

import nrpys

from antismash.common.html_renderer import Markup
from antismash.common.path import find_latest_database_version
from antismash.config import ConfigType
from antismash.detection.nrps_pks_domains import ModularDomain

from .data_structures import Prediction
from .name_mappings import get_substrate_by_name, SubstrateName
from .signatures import get_a_dom_signatures

SIGNATURE_FILE_NAME = "signatures.tsv"
MINIMUM_VERSION = "1.1"


def _get_signature_path(config: ConfigType, minimum_version: Optional[str] = None) -> str:
    """ A helper to construct the absolute path to the Stachelhaus signature file in the
        data directory.

        Arguments:
            config: the antiSMASH config

        Returns:
            the absolute path of the signature file
    """
    root = os.path.join(config.database_dir, "nrps_pks", "stachelhaus")
    version = find_latest_database_version(root)
    if minimum_version:
        min_major, min_minor = map(int, minimum_version.split("."))
        version_major, version_minor = map(int, version.split("."))
        if version_major < min_major or (version_major == min_major and version_minor < min_minor):
            return os.path.join(root, minimum_version, SIGNATURE_FILE_NAME)
    return os.path.join(root, version, SIGNATURE_FILE_NAME)


def _get_model_dir(config: ConfigType) -> str:
    """ A helper to construct the absolute path to the NRPS SVM model base dir in the
        data directory.

        Arguments:
            config: the antiSMASH config

        Returns:
            the absolute path of the NRPS SVM model base dir
    """
    root = os.path.join(config.database_dir, "nrps_pks", "svm")
    version = find_latest_database_version(root)
    return os.path.join(root, version)


def check_prereqs(options: ConfigType) -> list[str]:
    """ Check if the signatures file is present. """
    failure_messages: list[str] = []
    if "mounted_at_runtime" in options.database_dir:  # can't prepare this one
        return failure_messages
    datafile_path = _get_signature_path(options, minimum_version=MINIMUM_VERSION)
    if not os.path.exists(datafile_path):
        failure_messages.append(f"Failed to locate {datafile_path}")
    model_dir = _get_model_dir(options)
    if not os.path.exists(model_dir):
        failure_messages.append(f"Failed to locate {model_dir}")
    return failure_messages


def _get_norine_if_not_x(substrate: SubstrateName) -> str:
    if substrate.norine == "X":
        return substrate.short
    return substrate.norine


@dataclass
class StachelhausMatch:
    """A Stachelhaus-table-based A domain substrate prediction"""
    substrates: list[SubstrateName]
    signature: str
    aa10_score: float
    aa34_score: float

    def __str__(self) -> str:
        name = " or ".join(map(_get_norine_if_not_x, self.substrates))
        return f"{name} ({self.signature})"

    def html(self) -> str:
        """ Return HTML-friendly stringification """
        name = " or ".join(map(_get_norine_if_not_x, self.substrates))
        return (f'<dd>{name} <span class="serif">{self.signature}</span> '
                f'({round(self.aa34_score * 100)}% 8Å match)</dd>')

    def to_json(self) -> dict[str, Any]:
        """ Creates a JSON representation of a StachelhausMatch """
        return {
            "substrates": list(map(lambda x: x.to_json(), self.substrates)),
            "signature": self.signature,
            "aa10_score": self.aa10_score,
            "aa34_score": self.aa34_score,
        }

    @classmethod
    def from_json(cls, data: dict[str, Any]) -> "StachelhausMatch":
        """ Create a StachelhausMatch from the JSON representation """
        substrates = [SubstrateName.from_json(d) for d in data["substrates"]]
        return cls(substrates, data["signature"], data["aa10_score"], data["aa34_score"])

    @classmethod
    def from_nrpys(cls, pred: nrpys.StachPrediction) -> "StachelhausMatch":
        """ Create a StachelhausMatch from the nrpys.StachPrediction"""
        substrates = [get_substrate_by_name(name) for name in pred.name.split("|")]
        return cls(substrates, pred.aa10_sig, pred.aa10_score, pred.aa34_score)


NRPSPRED_TO_AS = {
    "bht": "R-ohTyr",
    "dhb": "2,3-dohBza",
    "horn": "ohOrn",
    "phg": "Pgl",
}


@dataclass
class SvmPrediction:
    """A Stachelhaus-table-based A domain substrate prediction"""
    name: str
    score: float
    substrates: list[SubstrateName]

    def __str__(self) -> str:
        if self.name == "N/A":
            return self.name
        substrates = ", ".join(map(_get_norine_if_not_x, self.substrates))
        return f"{self.name} ({substrates})"

    @classmethod
    def from_nrpys(cls, predictions: list[nrpys.Prediction]) -> "SvmPrediction":
        """ Load from an nrpys.Prediction, converting substrate strings to SubstrateNames """
        if not predictions:
            return cls("N/A", 0.0, [])

        pred = predictions[0]

        name, substrates = _resolve_name_substrates(pred)

        return cls(name, pred.score, substrates)

    def to_json(self) -> dict[str, Any]:
        """ Creates a JSON representation of an SvmPrediction """
        return {
            "name": self.name,
            "score": self.score,
            "substrates": list(map(lambda x: x.to_json(), self.substrates)),
        }

    @classmethod
    def from_json(cls, data: dict[str, Any]) -> "SvmPrediction":
        """ Create an SvmPrediction from the JSON representation """
        substrates = [SubstrateName.from_json(d) for d in data["substrates"]]
        return cls(data["name"], data["score"], substrates)


def _resolve_name_substrates(pred: nrpys.Prediction) -> tuple[str, list[SubstrateName]]:
    substrates: list[SubstrateName] = []

    if pred.name == "hydrophilic":
        substrates = list(map(get_substrate_by_name, [
            "Arg", "Asp", "Glu", "Asn", "Lys", "Gln",
            "Orn", "Aad",
        ]))
        return (pred.name, substrates)
    if pred.name == "hydrophobic-aliphatic":
        substrates = list(map(get_substrate_by_name, [
            "Ala", "Gly", "Val", "Leu", "Ile", "Abu",
            "Iva", "Ser", "Thr", "Hpg", "dHpg", "Cys",
            "Pro", "Pip",
        ]))
        return (pred.name, substrates)
    if pred.name == "hydrophobic-aromatic":
        substrates = list(map(get_substrate_by_name, [
            "Phe", "Tyr", "2,3-dohBza", "Pgl", "R-ohTyr",
        ]))
        return (pred.name, substrates)

    for part in pred.name.split(","):
        name = NRPSPRED_TO_AS.get(part, part)
        substrates.append(get_substrate_by_name(name))

    if pred.name == "asp,asn,glu,gln,aad":
        name = "Aliphatic chain with H-bond donor"
    elif pred.name == "cys":
        name = "Polar, uncharged (aliphatic with -SH)"
    elif pred.name == "dhb,sal":
        name = "Hydroxy-benzoic acid derivates"
    elif pred.name == "gly,ala,val,leu,ile,abu,iva":
        name = "Apolar, aliphatic"
    elif pred.name == "orn,lys,arg":
        name = "Long positively charged side chain"
    elif pred.name == "phe,trp,phg,tyr,bht":
        name = "Aromatic side chain"
    elif pred.name == "pro,pip":
        name = "Cyclic aliphatic chain (polar NH2 group)"
    elif pred.name == "ser,thr,dhpg,hpg":
        name = "Aliphatic chain or phenyl group with -OH"
    elif pred.name == "dhpg,hpg":
        name = "Polar, uncharged (hydroxy-phenyl)"
    elif pred.name == "gly,ala":
        name = "Tiny, hydrophilic, transition to aliphatic"
    elif pred.name == "orn,horn":
        name = "Orn and hydroxy-Orn specific"
    elif pred.name == "phe,trp":
        name = "Unpolar aromatic ring"
    elif pred.name == "tyr,bht":
        name = "Polar aromatic ring"
    elif pred.name == "val,leu,ile,abu,iva":
        name = "Aliphatic, branched hydrophobic"
    else:
        name = pred.name

    return (name, substrates)


class PredictorSVMResult(Prediction):
    """ Holds all the relevant results from nrpys for a domain """
    def __init__(self,
                 aa34: str,
                 aa10: str,
                 stachelhaus_matches: list[StachelhausMatch],
                 physicochemical_class: SvmPrediction,
                 large_cluster: SvmPrediction,
                 small_cluster: SvmPrediction,
                 single_amino: SvmPrediction,
                 ) -> None:
        super().__init__("nrpys")
        self.aa34 = aa34
        self.aa10 = aa10

        self.stachelhaus_matches = stachelhaus_matches
        self.stachelhaus_matches.sort(key=lambda x: (x.aa10_score, x.aa34_score), reverse=True)

        self.stachelhaus_quality = 0.0
        if stachelhaus_matches:
            self.stachelhaus_quality = max(map(lambda x: x.aa10_score, stachelhaus_matches))
        self.physicochemical_class = physicochemical_class
        self.large_cluster = large_cluster
        self.small_cluster = small_cluster
        self.single_amino = single_amino
        self.uncertain = aa34.count("-") >= 10

    def _get_classification(self) -> list[SubstrateName]:
        # comparing number of stach matches (n) to which category of SVM prediction
        # was made, and also to whether the SVM registered being outside of applicability domain
        # < = take stach, ^ = take SVM, & = take intersection of both, . = neither
        #    n   | single small/large/physico outside
        #    10      <            <              <
        #     9      &            &              <
        #     8      ^            &              <
        #  <= 7      ^            ^              .
        classification: list[SubstrateName] = []

        def get_best_stach_results() -> list[SubstrateName]:
            """Get the best Stachelhaus matches by aa10/aa34 score"""
            preds = [self.stachelhaus_matches[0]]
            for match in self.stachelhaus_matches[1:]:
                if match.aa34_score == preds[0].aa34_score:
                    preds.append(match)
            return list(
                itertools.chain.from_iterable(map(lambda x: x.substrates, preds)))

        if self.uncertain:
            if self.stachelhaus_quality >= 0.8:
                classification.extend(get_best_stach_results())
            return classification

        def stach_intersection_with_best_group() -> set[SubstrateName]:
            """ Finds the intersection of stach predictions with the tightest
                group of SVM predictions. If no SVM prediction, returns stach preds instead.
            """
            stach_preds: set[SubstrateName] = set(
                itertools.chain.from_iterable(
                    map(lambda x: x.substrates, self.stachelhaus_matches)))
            for group in [self.small_cluster, self.large_cluster, self.physicochemical_class]:
                if group.name == "N/A":
                    continue
                return stach_preds.intersection(set(group.substrates))
            return stach_preds

        stach_preds: set[SubstrateName] = set(
            itertools.chain.from_iterable(map(lambda x: x.substrates, self.stachelhaus_matches)))
        if self.stachelhaus_quality == 1.0:
            # return the best stachelhaus match(es)
            if self.stachelhaus_matches:
                classification.extend(get_best_stach_results())
        elif self.single_amino.name != "N/A":
            if self.stachelhaus_quality == 0.9:
                if self.single_amino.substrates[0] in stach_preds:
                    classification.extend(self.single_amino.substrates)
            else:
                classification.extend(self.single_amino.substrates)
        elif self.stachelhaus_quality >= 0.8:
            classification.extend(stach_intersection_with_best_group())
        else:  # < 8 and not uncertain
            for group in [self.single_amino, self.small_cluster,
                          self.large_cluster, self.physicochemical_class]:
                if group.name != "N/A":
                    classification.extend(group.substrates)
                    break
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
            note = "<strong>NOTE: uncertain match</strong><br>\n"
        qualifier = "weak"
        if self.stachelhaus_quality == 1.0:
            qualifier = "strong"
        elif self.stachelhaus_quality > 0.7:
            qualifier = "moderate"

        stachelhaus_results = "\n".join(map(lambda x: x.html(), self.stachelhaus_matches))

        raw = ("\n"
               "<dl><dt>SVM prediction details:</dt>\n"
               " <dd>"
               f"  {note}"
               "  <dl>"
               "   <dt>Predicted physicochemical class:</dt>\n"
               f"   <dd>{self.physicochemical_class}</dd>\n"
               "   <dt>Large clusters prediction:</dt>\n"
               f"   <dd>{self.large_cluster}</dd>\n"
               "   <dt>Small clusters prediction:</dt>\n"
               f"   <dd>{self.small_cluster}</dd>\n"
               "   <dt>Single AA prediction:</dt>\n"
               f"   <dd>{self.single_amino}</dd>\n"
               "  </dl>\n"
               " </dd>\n"
               "</dl>\n"
               "<dl><dt>Stachelhaus prediction details:</dt>\n"
               " <dd>\n"
               "  <dl>\n"
               "   <dt>Stachelhaus sequence:</dt>\n"
               f"   <dd><span class=\"serif\">{self.aa10}</span></dd>\n"
               "   <dt>Nearest Stachelhaus code(s):</dt>\n"
               f"   {stachelhaus_results}\n"
               "   <dt>Stachelhaus code match:</dt>\n"
               f"   <dd>{round(self.stachelhaus_quality * 100)}% ({qualifier})</dd>\n"
               "  </dl>\n"
               " </dd>\n"
               "</dl>\n")
        return Markup(raw)

    @classmethod
    def from_nrpys(cls, a_domain: nrpys.ADomain) -> "PredictorSVMResult":
        """ Generates a PredictorSVMResult from an nrpys.ADomain """
        aa34 = a_domain.aa34
        aa10 = a_domain.aa10

        stachelhaus_matches: list[StachelhausMatch] = []
        seen: set[str] = set()
        for pred in a_domain.get_stachelhaus():
            if pred.name in seen:
                continue
            seen.add(pred.name)
            stachelhaus_matches.append(StachelhausMatch.from_nrpys(pred))

        physiochemical_class = SvmPrediction.from_nrpys(
            a_domain.get_best(nrpys.PredictionCategory.ThreeClusterV2))
        large_custer = SvmPrediction.from_nrpys(
            a_domain.get_best(nrpys.PredictionCategory.LargeClusterV2))
        small_custer = SvmPrediction.from_nrpys(
            a_domain.get_best(nrpys.PredictionCategory.SmallClusterV2))
        single_amino = SvmPrediction.from_nrpys(
            a_domain.get_best(nrpys.PredictionCategory.SingleV2))

        return cls(aa34, aa10, stachelhaus_matches, physiochemical_class, large_custer,
                   small_custer, single_amino)

    def __str__(self) -> str:
        return "PredictorSVMResult: " + str(vars(self))

    def to_json(self) -> dict[str, Any]:
        stach_matches = list(map(lambda x: x.to_json(), self.stachelhaus_matches))
        return {
            "aa10": self.aa10,
            "aa34": self.aa34,
            "stachelhaus_matches": stach_matches,
            "physiochemical_class": self.physicochemical_class.to_json(),
            "large_cluster": self.large_cluster.to_json(),
            "small_cluster": self.small_cluster.to_json(),
            "single_amino": self.single_amino.to_json()
        }

    @classmethod
    def from_json(cls, json: dict[str, Any]) -> "PredictorSVMResult":
        stach_matches = [StachelhausMatch.from_json(d) for d in json["stachelhaus_matches"]]
        return PredictorSVMResult(json["aa10"], json["aa34"], stach_matches,
                                  SvmPrediction.from_json(json["physiochemical_class"]),
                                  SvmPrediction.from_json(json["large_cluster"]),
                                  SvmPrediction.from_json(json["small_cluster"]),
                                  SvmPrediction.from_json(json["single_amino"]),
                                  )


def run_nrpys(a_domains: list[ModularDomain], options: ConfigType) -> dict[str, Prediction]:
    """ Runs nrpys over the provided A domains.

        Arguments:
            a_domains: a list of ModularDomains, one for each A domain
            options: antismash options

        Returns:
            a dictionary mapping each domain name to a PredictorSVMResult
    """
    # NRPSPredictor: extract AMP-binding + 120 residues N-terminal of this domain,
    # extract 8 Angstrom residues and insert this into NRPSPredictor
    config = nrpys.Config()
    config.model_dir = _get_model_dir(options)
    config.stachelhaus_signatures = _get_signature_path(options)
    config.skip_v1 = True
    config.skip_v3 = True

    names: list[str] = []
    signatures: list[str] = []
    for sig, domain in [(get_a_dom_signatures(a_domain)[1], a_domain) for a_domain in a_domains]:
        if not sig:
            continue
        names.append(domain.get_name())
        signatures.append(sig)

    results: dict[str, Prediction] = {}
    for res in nrpys.run(config, names, signatures):
        results[res.name] = PredictorSVMResult.from_nrpys(res)

    return results
