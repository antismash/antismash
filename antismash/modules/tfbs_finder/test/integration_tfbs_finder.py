# License: GNU Affero General Public License v3 or later
# A copy of GNU AGPL v3 should have been included in this software package in LICENSE.txt.

# for test files, silence irrelevant and noisy pylint warnings
# pylint: disable=protected-access,missing-docstring

import unittest

from Bio.Seq import Seq

import antismash
from antismash.common.test import helpers
from antismash.config import destroy_config, build_config
from antismash.modules import tfbs_finder


class TestTFBSFinder(unittest.TestCase):
    def setUp(self):
        self.options = build_config(["--minimal", "--enable-html", "--tfbs"],
                                    isolated=True, modules=antismash.get_all_modules())

    def tearDown(self):
        destroy_config()

    def test_options(self):
        assert not tfbs_finder.check_options(self.options)
        assert not tfbs_finder.check_prereqs(self.options)
        assert tfbs_finder.is_enabled(self.options)

        options = build_config(["--minimal", "--enable-html", "--tfbs",
                                "--tfbs-pvalue", "-1", "--tfbs-range=-10:-10"],
                               isolated=True, modules=antismash.get_all_modules())
        issues = tfbs_finder.check_options(options)
        assert len(issues) == 2
        assert "TFBS finder p-value is negative" in issues[0]
        assert "TFBS finder range is negative" in issues[1]

    def test_balhimycin(self):
        genbank = helpers.get_path_to_balhymicin_genbank()
        results = helpers.run_and_regenerate_results_for_module(genbank, tfbs_finder, self.options)
        assert isinstance(results, tfbs_finder.TFBSFinderResults)
        assert results.pvalue == self.options.tfbs_pvalue
        assert (results.upstream_range, results.start_overlap) == \
            tfbs_finder.tfbs_finder.parse_range(self.options.tfbs_range)
        assert results.record_id == 'Y16952.3'
        assert len(results.hits_by_region) == 1

        hits = results.get_hits_by_region(1)
        assert len(hits) == 17
        strong = results.get_hits_by_region(1, confidence=tfbs_finder.tfbs_finder.Confidence.STRONG)
        assert len(strong) == 9
        assert not results.get_hits_by_region(1, confidence=tfbs_finder.tfbs_finder.Confidence.WEAK)

        best = max(hits, key=lambda hit: hit.score)
        assert best.name == "VqsM"
        assert best.start == 1398
        assert best.strand == 1
        self.assertAlmostEqual(best.score, 20.469, delta=1e-3)

    def test_run_moods(self):
        sequence = Seq("GAATAACGTTAGGCTCAACGTTGCT")
        bg = tfbs_finder.tfbs_finder.get_bg_distribution(sequence)
        matrices = tfbs_finder.tfbs_finder.load_matrices(tfbs_finder.tfbs_finder.PWM_PATH)
        region_start = 100
        all_hits = tfbs_finder.tfbs_finder.run_moods(sequence, bg, matrices, 0.01, region_start)
        matrix_hits = {matrix.name: hits for matrix, hits in zip(matrices, all_hits) if hits}
        assert len(matrix_hits) == 45
        # the scores will be largely useless thanks to the high threshold, but
        # the positional information is worth testing
        expected_positions = {
            "BldD": [(10, -1)],
            "CsoR": [(2, -1)],
            "DmdR1": [(3, 1)],
            "OsdR": [(10, 1)],
        }
        for name, hits in expected_positions.items():
            for hit, expected in zip(matrix_hits[name], hits):
                assert (hit.pos - region_start, hit.strand) == expected
