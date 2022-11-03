# License: GNU Affero General Public License v3 or later
# A copy of GNU AGPL v3 should have been included in this software package in LICENSE.txt.

# for test files, silence irrelevant and noisy pylint warnings
# pylint: disable=use-implicit-booleaness-not-comparison,protected-access,missing-docstring

import unittest

from antismash.common import path, fasta
from antismash.config import build_config, destroy_config
from antismash.modules.nrps_pks.kr_analysis.kr_analysis import run_kr_analysis


class TestKRAnalysis(unittest.TestCase):
    def setUp(self):
        build_config([])
        self.query_data = fasta.read_fasta(path.get_full_path(__file__, 'data', 'SCO_genes.fasta'))

    def tearDown(self):
        destroy_config()

    def test_full_run(self):
        active, stereo = run_kr_analysis(self.query_data)
        active = {key: val.prediction for key, val in active.items()}
        assert active == {'SCO0111_KR1': "inactive",
                          'SCO0126_KR1': "active",
                          'SCO0256_KR1': "inactive",
                          'SCO5086_KR1': "inactive",
                          'SCO5097_KR1': "inactive",
                          'SCO6249_KR1': "inactive",
                          'SCO6264_KR1': "inactive",
                          'SCO6273_KR1': "active",
                          'SCO6274_KR1': "active",
                          'SCO6274_KR2': "active",
                          'SCO6275_KR1': "active",
                          'SCO6275_KR2': "active",
                          'SCO6282_KR1': "inactive",
                          'SCO6827_KR1': "active",
                          'SCO6829_KR1': "inactive"}
        stereo = {key: val.prediction for key, val in stereo.items()}
        assert stereo == {'SCO0111_KR1': 'C2',
                          'SCO0256_KR1': 'C1',
                          'SCO5086_KR1': 'C1',
                          'SCO5097_KR1': 'C2',
                          'SCO6249_KR1': 'C2',
                          'SCO6264_KR1': 'C2',
                          'SCO6273_KR1': 'B1',
                          'SCO6274_KR1': 'B1',
                          'SCO6275_KR1': 'B1',
                          'SCO6275_KR2': 'B1',
                          'SCO6282_KR1': 'C2',
                          'SCO6829_KR1': 'C2'}
