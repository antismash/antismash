# License: GNU Affero General Public License v3 or later
# A copy of GNU AGPL v3 should have been included in this software package in LICENSE.txt.

"""
Standalone Stachelhaus code lookup-table-based A domain substrate predictions
"""

from collections import defaultdict
from dataclasses import dataclass
from typing import Any, Optional
import os

from antismash.common.html_renderer import Markup
from antismash.common.path import find_latest_database_version
from antismash.config import ConfigType
from antismash.detection.nrps_pks_domains import ModularDomain

from .data_structures import Prediction
from .name_mappings import SubstrateName, get_substrate_by_name
from .signatures import get_a_dom_signatures

SIGNATURE_FILE_NAME = "signatures.tsv"


def check_prereqs(options: ConfigType) -> list[str]:
    """ Check if the signatures file is present. """
    failure_messages: list[str] = []
    if "mounted_at_runtime" in options.database_dir:  # can't prepare this one
        return failure_messages
    datafile_path = _get_datafile_path(options)
    if not os.path.exists(datafile_path):
        failure_messages.append(f"Failed to locate {datafile_path}")
    return failure_messages


@dataclass
class StachSignature:
    """ A Stachelhaus signature hit """
    aa10: str
    aa34: str
    substrates: set[SubstrateName]
    smiles_key: str

    def to_json(self) -> dict[str, Any]:
        """ Creates a JSON representation of a StachSignature """
        data = {
            "aa10": self.aa10,
            "aa34": self.aa34,
            "substrates": sorted(map(lambda x: x.to_json(), self.substrates),
                                 key=lambda x: x["short"]),
            "smiles_key": self.smiles_key
        }
        return data

    @classmethod
    def from_json(cls, data: dict[str, Any]) -> "StachSignature":
        """ Creates a StachSignature from a JSON representation """
        substrates = set(map(SubstrateName.from_json, data["substrates"]))
        return cls(data["aa10"], data["aa34"], substrates, data["smiles_key"])


def _get_datafile_path(config: ConfigType) -> str:
    """ A helper to construct absolute paths to files in the knownclusterblast
        data directory.

        Arguments:
            filename: the name only of the file

        Returns:
            the absolute path of the file
    """
    root = os.path.join(config.database_dir, "nrps_pks", "stachelhaus")
    version = find_latest_database_version(root)
    return os.path.join(root, version, SIGNATURE_FILE_NAME)


KNOWN_STACH_CODES: dict[str, list[StachSignature]] = {}


def init_data(config: ConfigType) -> dict[str, list[StachSignature]]:
    """ Builds a mapping of Stachelhaus prediction to code from the PARAS data """
    if KNOWN_STACH_CODES:
        return KNOWN_STACH_CODES

    data_file = _get_datafile_path(config)

    mappings: dict[str, list[StachSignature]] = defaultdict(list)

    with open(data_file, encoding="utf-8") as handle:
        for line in handle:
            parts = line.strip().split('\t')
            assert len(parts) == 5, f"malformed Stachelhaus signature file {parts}"
            aa10, aa34, raw_names, smiles_key, _ = parts
            assert len(aa10) == 10, "malformed Stachelhaus signature file"
            assert "-" not in aa10, "malformed Stachelhaus signature"
            names = [get_substrate_by_name(name) for name in raw_names.split("|")]
            assert len(names) > smiles_key.count("|"), "malformed Stachelhaus signature file"
            mappings[aa10].append(StachSignature(aa10, aa34, set(names), smiles_key))

    KNOWN_STACH_CODES.update(mappings)
    return mappings


class StachelhausPrediction(Prediction):
    """ Holds all the relevant info for a Stachelhaus signature prediction for a domain """
    def __init__(self, signature: Optional[StachSignature] = None,
                 original_signature: Optional[str] = None) -> None:
        super().__init__("Stachelhaus")
        self.signature = signature
        self.aa10 = original_signature if not signature else signature.aa10
        if not self.aa10:
            raise ValueError("Need to provide either a StachSignature or an original 10 AA signature")

    def get_classification(self, as_norine: bool = False) -> list[str]:
        classification = []
        if self.signature:
            for substrate in self.signature.smiles_key.split("|"):
                sub = get_substrate_by_name(substrate)
                if as_norine:
                    classification.append(sub.norine)
                else:
                    classification.append(sub.short)
        return classification

    # TODO: Also show the aa10 sequence and all the substrate names
    def as_html(self) -> Markup:
        header = "Stachelhaus lookup:"
        prediction = "No exact match."
        footer = f"10 AA signature {self.aa10}"
        if self.signature:
            names = []
            for substrate in self.signature.smiles_key.split("|"):
                names.append(get_substrate_by_name(substrate).long)
            prediction = " or ".join(names)
        return Markup(f"{header} {prediction}<br>{footer}")

    def to_json(self) -> dict[str, Any]:
        signature = None
        if self.signature:
            signature = self.signature.to_json()
        return {"signature": signature, "aa10": self.aa10}

    @classmethod
    def from_json(cls, json: dict[str, Any]) -> "Prediction":
        if json["signature"]:
            return cls(StachSignature.from_json(json["signature"]))
        return cls(original_signature=json["aa10"])


def run_stachelhaus(a_domains: list[ModularDomain], config: ConfigType) -> dict[str, Prediction]:
    """ Compares the extracted Stachelhaus sequences against the table of known sequences """
    predictions: dict[str, Prediction] = {}
    signatures = [get_a_dom_signatures(a_domain) for a_domain in a_domains]
    mappings = init_data(config)

    for (aa10, aa34), a_domain in zip(signatures, a_domains):
        # if we didn't get signatures, it's not really an A domain, don't predict anything
        if not (aa10 and aa34):
            continue

        options = mappings.get(aa10, [])
        # if there's no match, record it, check the next domain
        if not options:
            predictions[a_domain.get_name()] = StachelhausPrediction(original_signature=aa10)
            continue

        # if there's exactly one match, record it, check the next domain
        if len(options) == 1:
            predictions[a_domain.get_name()] = StachelhausPrediction(options[0])
            continue

        consensus_signature = pick_winning_substrate(aa10, aa34, options)

        predictions[a_domain.get_name()] = StachelhausPrediction(consensus_signature)
    return predictions


def pick_winning_substrate(aa10: str, aa34: str, options: list[StachSignature]) -> StachSignature:
    """ Pick the substrate(s) that occur(s) most often as the winner(s) """
    # check if there's an exact 34 AA signature match
    # but also get the smiles_key counts to build a consensus
    # and keep track of substrates in general
    aa34_match = False
    counts: dict[str, int] = defaultdict(int)
    substrates: dict[str, set[SubstrateName]] = defaultdict(set)
    option: Optional[StachSignature] = None
    for option in options:
        counts[option.smiles_key] += 1
        substrates[option.smiles_key].update(option.substrates)
        if aa34 == option.aa34:
            aa34_match = True
            break

    # if there's an exact 34 AA signature match, record it, check the next domain
    if aa34_match:
        return option

    keys_by_count = sorted(counts.items(), key=lambda x: x[1], reverse=True)
    consensus = keys_by_count[0][0]
    consensus_substrates = substrates[consensus]
    for key, count in keys_by_count[1:]:
        if count < keys_by_count[0][1]:
            break
        consensus += f"|{key}"
        consensus_substrates.update(substrates[key])

    return StachSignature(aa10, aa34, consensus_substrates, consensus)
