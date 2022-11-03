# License: GNU Affero General Public License v3 or later
# A copy of GNU AGPL v3 should have been included in this software package in LICENSE.txt.

# for test files, silence irrelevant and noisy pylint warnings
# pylint: disable=use-implicit-booleaness-not-comparison,protected-access,missing-docstring

import unittest

from antismash.modules.thiopeptides import rodeo


class TestRodeo(unittest.TestCase):
    def test_acquire_rodeo_heuristics(self):
        """Test thiopeptides.acquire_rodeo_heuristics()"""
        leader = "MSDITASRVESLDLQDLDLSELTVTSLRDTVALPENGA"
        core = "SWGSCSCQASSSCAQPQDM"
        domains = [
            'LANC_like',
            'Pkinase',
            'Lant_dehyd_N',
            'Lant_dehyd_C',
            'Lant_dehydr_C',
            'YcaO',
            'PF00881',
            'Lant_dehydr_C'
        ]
        expected_score = 15
        expected_tabs = [1, 1, 1, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 1, 0, 0,
                         0, 0, 0, 0, 0, 0, 0, 1, 1, 0, 1, 1, 0, 0, 0, 0, 1]
        print(expected_tabs)
        score, tabs = rodeo.acquire_rodeo_heuristics(leader, core, domains)
        print(tabs)
        self.assertEqual(expected_score, score)
        self.assertEqual(expected_tabs, tabs)

        core = "SWGSCSCQASSSCAQPQDMX"
        score, tabs = rodeo.acquire_rodeo_heuristics(leader, core, domains)
        self.assertEqual(expected_score, score)
        self.assertEqual(expected_tabs, tabs)

    def test_run_rodeo_svm(self):
        # from the thiopeptide.gbk test file, BGC0001155.1
        inputs = [1, 0, 1, 1, 1, 1, 1, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0,
                  0, 0, 0, 0, 0, 1, 1, 0, 1, 1, 0, 0, 0, 0, 1, 1, 2, 0, 0, 2,
                  3.5, 5734.265325999999, 4182.791014999999, 1569.4848759999998,
                  55, 39, 16, 0.41025641025641024, 0.625, 3, 0, 3, 2, 0, 0, 6,
                  4, 1, 0, 4, 0, 6, 1, 2, 3, 1, 0, 0, 3, 1, 9, 0, 9, 20, 4, 1,
                  0, 0, 1, 6, 0, 0, 1, 0, 0, 0, 0, 0, 1, 1, 4, 0, 0, 0, 1, 1,
                  0, 0, 0, 3, 4, 4, 0, 3, 3, 6, 0, 6, 5, 1, 0, 4, 0, 6, 2, 3,
                  7, 1, 0, 0, 4, 2, 9, 0, 9, 23, 8]
        res = rodeo.run_rodeo_svm(inputs)
        assert res == 10

        # from CP009500.1, a thiopeptide cluster
        inputs = [1, 0, 1, 1, 1, 1, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
                  0, 0, 0, 1, 0, 0, 0, 0, 0, 1, 0, 0, 1, 1, 0, 0, 0, 0, 0, 0,
                  0.0, 15234.355701999995, 7588.556352999996, 7663.809914000001,
                  130, 64, 66, 1.03125, 0.09090909090909091, 3, 5, 6, 0, 0, 3,
                  4, 6, 2, 1, 4, 0, 2, 6, 4, 3, 4, 3, 4, 4, 13, 10, 5, 15, 20,
                  7, 10, 8, 6, 2, 0, 2, 4, 1, 4, 2, 9, 0, 4, 1, 3, 4, 2, 2, 1,
                  1, 4, 10, 8, 18, 27, 6, 13, 13, 12, 2, 0, 5, 8, 7, 6, 3, 13,
                  0, 6, 7, 7, 7, 6, 5, 5, 5, 17, 20, 13, 33, 47, 13]
        res = rodeo.run_rodeo_svm(inputs)
        assert res == 0

        # from CP003082.1, a thiopeptide cluster
        inputs = [1, 0, 1, 1, 0, 1, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
                  1, 0, 0, 1, 0, 0, 0, 0, 0, 1, 1, 0, 1, 1, 0, 2, 0, 0, 4, 2,
                  2.0, 20322.627290999986, 3456.9570819999994, 16883.680774000008,
                  192, 33, 159, 4.818181818181818, 0.11320754716981132, 4, 1, 0,
                  1, 0, 0, 0, 3, 0, 1, 5, 0, 1, 0, 0, 2, 3, 2, 1, 9, 3, 0, 1, 1,
                  23, 5, 16, 6, 16, 8, 0, 7, 5, 14, 2, 13, 15, 1, 2, 5, 5, 8,
                  10, 2, 3, 21, 10, 21, 7, 28, 81, 18, 20, 7, 16, 9, 0, 7, 5,
                  17, 2, 14, 20, 1, 3, 5, 5, 10, 13, 4, 4, 30, 13, 21, 8, 29,
                  104, 23]
        res = rodeo.run_rodeo_svm(inputs)
        assert res == 0

    def test_statistics(self):
        stats = rodeo.ThioStatistics("CCSSSTTTT")
        assert stats.c_repeats == 2
        assert stats.s_repeats == 3
        assert stats.t_repeats == 4
        assert stats.cst_repeats == 9
        assert stats.block_repeats == 3
        assert stats.heteroblocks == 1
        assert stats.average_heteroblock_length == 9

        stats = rodeo.ThioStatistics("CCSxSSTxTTT")
        assert stats.c_repeats == 2
        assert stats.s_repeats == 2
        assert stats.t_repeats == 3
        assert stats.cst_repeats == 7
        assert stats.block_repeats == 3
        assert stats.heteroblocks == 3
        assert stats.average_heteroblock_length == 3

    def test_run_rodeo(self):
        # from CP003082
        leader = "MLTLAAVGRSAVSVVIWVWVGLAVLTYLGVNVT"
        core = ("PVVASAGVVGVALGFGAQSLVKDFLSGIFMLIEDQYGVGDTIDVGDVVGTVEDVSLRLTTLR"
                "DIHGTQWFVRNGEILRIGNFSQDYAVALINIPIALDENASDAINAVTDAVNAAAQEPAINDV"
                "LLDSPTVDGVNSIGLDHMLIRARVTTLPDQQWYVT")
        domains = {'Peptidase_S8', 'PF00881', 'YcaO', 'Lant_dehydr_C', 'PF00005', 'Lant_dehyd_N'}
        valid, score = rodeo.run_rodeo(leader, core, domains)
        assert not valid
        assert score == 6

        # from BGC0001155.1
        leader = "MSEMELNLNDLPMDVFEMADSGMEVESLTAGHGMPEVGA"
        core = "SCNCVCGFCCSCSPSA"
        domains = {'PF04055', 'p450', 'Lant_dehyd_N', 'thio_amide',
                   'Lant_dehydr_C', 'PF00881', 'Lant_dehyd_C', 'YcaO'}
        valid, score = rodeo.run_rodeo(leader, core, domains)
        assert valid
        assert score == 27
