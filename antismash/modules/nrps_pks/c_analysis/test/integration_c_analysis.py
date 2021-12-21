# License: GNU Affero General Public License v3 or later
# A copy of GNU AGPL v3 should have been included in this software package in LICENSE.txt.

# for test files, silence irrelevant and noisy pylint warnings
# pylint: disable=protected-access,missing-docstring

import unittest

from antismash.common import path, fasta
from antismash.config import build_config, destroy_config
from antismash.modules.nrps_pks.c_analysis.c_analysis import run_c_analysis


class TestCAnalysis(unittest.TestCase):
    def setUp(self):
        build_config([])
        c_data = fasta.read_fasta(path.get_full_path(__file__, 'data', 'C_domains.fasta'))
        self.c_data = [(name.split()[0], name.split()[1], seq)
                            for (name, seq) in c_data.items()]

    def tearDown(self):
        destroy_config()

    def test_full_run(self):
        activity, activesite = run_c_analysis(self.c_data).values()
        activity = {key: val.prediction for key, val in activity.items()}
        assert activity == { # Vancomycin
                          'vcm2_C1': "active",
                          'vcm2_C2': "active",
                          'vcm2_E1': "active",
                          'vcm4_X1': "inactive",  # Inactive as a C domain
                          # Vibriobactin
                          'vabF_C2': "active",
                          'vabF_C3': "inactive",
                          'vabF_C4': "active",
                          # Anguibactin
                          'angR_C1': "inactive",
                          # Cupriachelin
                          'cucG_C1': "inactive",
                          # Histicorrugatin
                          'hcsF_C1': "inactive",
                          'hcsH_C2': "inactive",
                          # Arthrofactin
                          'arfA_C1': "active",
                          'arfA_C2': "active"}
        activesite = {key: val.prediction for key, val in activesite.items()}
        assert activesite == {'vcm2_C1': "HHILVDG",
                              'vcm2_C2': "HHILLDG",
                              'vcm2_E1': "HHLSVDG",
                              'vcm4_X1': "HRILADD",
                              'vabF_C2': "DALIVDG",
                              'vabF_C3': "HQIISEQ",
                              'vabF_C4': "HHIVLDE",
                              'angR_C1': "NSVVVDN",
                              'cucG_C1': "SPLAADR",
                              'hcsF_C1': "ASSHFDL",
                              'hcsH_C2': "HALRLDT",
                              'arfA_C1': "HHLIVDG",
                              'arfA_C2': "HHLAMDH"}
