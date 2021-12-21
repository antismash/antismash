# License: GNU Affero General Public License v3 or later
# A copy of GNU AGPL v3 should have been included in this software package in LICENSE.txt.

# for test files, silence irrelevant and noisy pylint warnings
# pylint: disable=protected-access,missing-docstring

import unittest
from unittest.mock import patch, MagicMock

from antismash.config import build_config, destroy_config
from antismash.modules.nrps_pks.c_analysis.c_analysis import is_active, run_c_analysis


class TestCAnalysis(unittest.TestCase):
    def setUp(self):
        build_config([])
        self.patched_gca = patch("antismash.common.brawn.get_cached_alignment")
        self.patched_gap = patch("antismash.common.brawn.get_aligned_pair")
        self.patched_ebrp = patch("antismash.common.utils.extract_by_reference_positions")

    def tearDown(self):
        destroy_config()
        patch.stopall()

    def test_is_active(self):
        data = [
            # "regular" C domains
            ("HHILVDG", "Cglyc", "active"),
            ("HHLSVDG", "Epimerization", "active"),
            ("HRILADD", "X", "inactive"),
            # Heterocyclization special case
            ("DALIVGG", "Heterocyclization", "active"),
            ("HQIISEQ", "Heterocyclization", "inactive"),
            # Condensation_sid special case
            ("HLLLWDG", "Condensation_sid", "active"),
            ("MAGICXX", "Condensation_sid", "inactive"),
        ]

        for signature, c_type, expected in data:
            result = is_active(signature, c_type)
            assert result == expected, f"{signature}|{c_type}: expected {expected!r}, got {result!r}"

    def test_run_c_analysis(self):
        signature = "HHILVDG"
        self.patched_gca.start()
        mock_pair = self.patched_gap.start()
        mock_pair.return_value = (MagicMock, MagicMock)
        mock_extract = self.patched_ebrp.start()
        mock_extract.side_effect = [signature, None]
        fake_data = [
            ("working", "Condensation", "MAGICHAT"),
            ("broken", "Condensation", "MAGICCAT"),
        ]
        results = run_c_analysis(fake_data)

        def get_classification(tool , name):
            return results[tool][name].get_classification()[0]


        assert get_classification('c_activity', 'working') == "active"
        assert get_classification('c_activesite', 'working') == signature
        assert get_classification('c_activity', 'broken') == "unknown"

        assert results['c_activesite']
