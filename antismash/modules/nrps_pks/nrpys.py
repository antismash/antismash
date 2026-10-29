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

from antismash.common.html_renderer import Markup, collapser_end, collapser_start
from antismash.common.path import find_latest_database_version
from antismash.config import ConfigType
from antismash.detection.nrps_pks_domains import ModularDomain

from .data_structures import Prediction
from .name_mappings import get_substrate_by_name, Substrate, Category, lookup
from .signatures import get_a_dom_signatures

SIGNATURE_FILE_NAME = "signatures.tsv"
MINIMUM_STACH_VERSION = "1.1"
MINIMUM_SVM_VERSION = "3.0"


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


def _get_model_dir(config: ConfigType, minimum_version: str | None = None) -> str:
    """ A helper to construct the absolute path to the NRPS SVM model base dir in the
        data directory.

        Arguments:
            config: the antiSMASH config
            minimum_version: the minimum version of the NRPS SVM model to use

        Returns:
            the absolute path of the NRPS SVM model base dir
    """
    root = os.path.join(config.database_dir, "nrps_pks", "svm")
    version = find_latest_database_version(root)
    if minimum_version:
        min_major, min_minor = map(int, minimum_version.split("."))
        version_major, version_minor = map(int, version.split("."))
        if version_major < min_major or (version_major == min_major and version_minor < min_minor):
            return os.path.join(root, minimum_version)
    return os.path.join(root, version)


def check_prereqs(options: ConfigType) -> list[str]:
    """ Check if the signatures file is present. """
    failure_messages: list[str] = []
    if "mounted_at_runtime" in options.database_dir:  # can't prepare this one
        return failure_messages
    try:
        datafile_path = _get_signature_path(options, minimum_version=MINIMUM_STACH_VERSION)
        if not os.path.exists(datafile_path):
            failure_messages.append(f"Signature file not found '{datafile_path}'")
    except ValueError as err:
        failure_messages.append(f"nrpys: {err}")

    try:
        model_dir = _get_model_dir(options, minimum_version=MINIMUM_SVM_VERSION)
        if not os.path.exists(model_dir):
            failure_messages.append(f"Model directory not found '{model_dir}'")
    except ValueError as err:
        failure_messages.append(f"nrpys: {err}")

    return failure_messages


@dataclass
class StachelhausMatch:
    """A Stachelhaus-table-based A domain substrate prediction"""
    substrates: list[Substrate]
    signature: str
    aa10_score: float
    aa34_score: float

    def __str__(self) -> str:
        name = " or ".join([substrate.short for substrate in self.substrates])
        return f"{name} ({self.signature})"

    def html(self) -> str:
        """ Return HTML-friendly stringification """
        name = " or ".join([substrate.short for substrate in self.substrates])
        return (f'<dd>{name} <span class="serif">{self.signature}</span> '
                f'({round(self.aa34_score * 100)}% 8Å match)</dd>')

    def to_json(self) -> dict[str, Any]:
        """ Creates a JSON representation of a StachelhausMatch """
        return {
            "substrates": [x.to_json() for x in self.substrates],
            "signature": self.signature,
            "aa10_score": self.aa10_score,
            "aa34_score": self.aa34_score,
        }

    @classmethod
    def from_json(cls, data: dict[str, Any]) -> "StachelhausMatch":
        """ Create a StachelhausMatch from the JSON representation """
        substrates = [Substrate.from_json(d) for d in data["substrates"]]
        return cls(substrates, data["signature"], data["aa10_score"], data["aa34_score"])

    @classmethod
    def from_nrpys(cls, pred: nrpys.StachPrediction) -> "StachelhausMatch":
        """ Create a StachelhausMatch from the nrpys.StachPrediction"""
        substrates = [get_substrate_by_name(name) for name in pred.name.split("|")]
        return cls(substrates, pred.aa10_sig, pred.aa10_score, pred.aa34_score)


@dataclass
class SvmPrediction:
    """A Stachelhaus-table-based A domain substrate prediction"""
    name: str
    score: float
    substrates: list[Substrate]

    def __str__(self) -> str:
        if self.name == "N/A":
            return self.name
        substrates = ", ".join([substrate.short for substrate in self.substrates])
        return f"{self.name}: {substrates}"

    def as_html(self) -> Markup:
        """ Return HTML-friendly stringification """
        if self.name == "N/A":
            return Markup(self.name)
        substrates = ", ".join([substrate.short for substrate in self.substrates])
        return Markup(f"{self.name}") + \
                collapser_start("svm-pred", level="none") + \
                Markup(f"{substrates}") +\
                collapser_end()

    @classmethod
    def from_nrpys(cls, predictions: list[nrpys.Prediction]) -> "SvmPrediction":
        """ Load from an nrpys.Prediction, converting substrate strings corresponding objects """
        if not predictions:
            return cls("N/A", 0.0, [])

        pred = predictions[0]

        # predictions can be a category or a substrate
        category_or_substrate = lookup(pred.name)
        if isinstance(category_or_substrate, Category):
            name = category_or_substrate.name
            substrates = list(category_or_substrate.substrates)
        else:
            name = category_or_substrate.short
            substrates = [category_or_substrate]

        return cls(name, pred.score, substrates)

    def to_json(self) -> dict[str, Any]:
        """ Creates a JSON representation of an SvmPrediction """
        return {
            "name": self.name,
            "score": self.score,
            "substrates": [x.to_json() for x in self.substrates],
        }

    @classmethod
    def from_json(cls, data: dict[str, Any]) -> "SvmPrediction":
        """ Create an SvmPrediction from the JSON representation """
        substrates = [Substrate.from_json(d) for d in data["substrates"]]
        return cls(data["name"], data["score"], substrates)


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

    def _get_classification(self) -> list[Substrate]:
        # comparing number of stach matches (n) to which category of SVM prediction
        # was made, and also to whether the SVM registered being outside of applicability domain
        # < = take stach, ^ = take SVM, & = take intersection of both, . = neither
        #    n   | single small/large/physico outside
        #    10      <            <              <
        #     9      &            &              <
        #     8      ^            &              <
        #  <= 7      ^            ^              .
        classification: list[Substrate] = []

        def get_best_stach_results() -> list[Substrate]:
            """Get the best Stachelhaus matches by aa10/aa34 score"""
            preds = [self.stachelhaus_matches[0]]
            for match in self.stachelhaus_matches[1:]:
                if match.aa34_score == preds[0].aa34_score:
                    preds.append(match)
            return list(
                itertools.chain.from_iterable(x.substrates for x in preds))

        if self.uncertain:
            if self.stachelhaus_quality >= 0.8:
                classification.extend(get_best_stach_results())
            return classification

        def stach_intersection_with_best_group() -> set[Substrate]:
            """ Finds the intersection of stach predictions with the tightest
                group of SVM predictions. If no SVM prediction, returns stach preds instead.
            """
            stach_preds: set[Substrate] = set(
                itertools.chain.from_iterable(
                    x.substrates for x in self.stachelhaus_matches))
            for group in [self.small_cluster, self.large_cluster, self.physicochemical_class]:
                if group.name == "N/A":
                    continue
                return stach_preds.intersection(set(group.substrates))
            return stach_preds

        stach_preds: set[Substrate] = set(
            itertools.chain.from_iterable(x.substrates for x in self.stachelhaus_matches))
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

        def mapper(substrate: Substrate) -> str:
            if as_norine:
                return substrate.norine
            return substrate.short

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

        stachelhaus_results = "\n".join(x.html() for x in self.stachelhaus_matches)

        raw = ("\n"
               "<dl><dt>SVM prediction details:</dt>\n"
               " <dd>"
               f"  {note}"
               "  <dl>"
               "   <dt>Predicted physicochemical class:</dt>\n"
               f"   <dd>{self.physicochemical_class.as_html()}</dd>\n"
               "   <dt>Large clusters prediction:</dt>\n"
               f"   <dd>{self.large_cluster.as_html()}</dd>\n"
               "   <dt>Small clusters prediction:</dt>\n"
               f"   <dd>{self.small_cluster.as_html()}</dd>\n"
               "   <dt>Single AA prediction:</dt>\n"
               f"   <dd>{self.single_amino.name}</dd>\n"
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
        assert len(aa10) == 10, aa10
        assert len(aa34) > 10, aa34

        stachelhaus_matches: list[StachelhausMatch] = []
        seen: set[str] = set()
        for pred in a_domain.get_stachelhaus():
            if pred.name in seen:
                continue
            seen.add(pred.name)
            stachelhaus_matches.append(StachelhausMatch.from_nrpys(pred))

        physiochemical_class = SvmPrediction.from_nrpys(
            a_domain.get_best(nrpys.PredictionCategory.ThreeClusterV3))
        large_custer = SvmPrediction.from_nrpys(
            a_domain.get_best(nrpys.PredictionCategory.LargeClusterV3))
        small_custer = SvmPrediction.from_nrpys(
            a_domain.get_best(nrpys.PredictionCategory.SmallClusterV3))
        single_amino = SvmPrediction.from_nrpys(
            a_domain.get_best(nrpys.PredictionCategory.SingleV3))

        return cls(aa34, aa10, stachelhaus_matches, physiochemical_class, large_custer,
                   small_custer, single_amino)

    def __str__(self) -> str:
        return "PredictorSVMResult: " + str(vars(self))

    def to_json(self) -> dict[str, Any]:
        stach_matches = [x.to_json() for x in self.stachelhaus_matches]
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
    config.skip_v2 = True

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
