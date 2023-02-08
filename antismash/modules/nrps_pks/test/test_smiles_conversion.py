# License: GNU Affero General Public License v3 or later
# A copy of GNU AGPL v3 should have been included in this software package in LICENSE.txt.

# for test files, silence irrelevant and noisy pylint warnings
# pylint: disable=use-implicit-booleaness-not-comparison,protected-access,missing-docstring

import os
import unittest

from antismash.common import path
from antismash.modules.nrps_pks.smiles_generator import get_all_smiles
from antismash.modules.nrps_pks.nrps_predictor import map_nrpspredicor_to_norine
from antismash.modules.nrps_pks.parsers import ALLOWABLE_PREDICTION_CHARACTERS

NRPS_PKS_DIR = os.path.dirname(__file__)
SIGNATURE_FILE = path.get_full_path(
    NRPS_PKS_DIR, "external", "NRPSPredictor2", "data", "labeled_sigs"
)


def load_sigs() -> set[str]:
    signatures: set[str] = set()
    with open(SIGNATURE_FILE, "r", encoding="utf-8") as handle:
        for line in handle:
            raw_name, _ = line.split('\t', 1)
            names = raw_name.split("/")
            for name in names:
                assert set(name).issubset(ALLOWABLE_PREDICTION_CHARACTERS)
                signatures.add(name)
    return signatures


SIGNATURES = load_sigs()


class TestSmilesConversion(unittest.TestCase):
    def test_all_signatures_have_smiles(self) -> None:
        smiles = get_all_smiles()
        for signature in SIGNATURES:
            assert map_nrpspredicor_to_norine(signature).lower() in smiles
