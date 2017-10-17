# License: GNU Affero General Public License v3 or later
# A copy of GNU AGPL v3 should have been included in this software package in LICENSE.txt.

# for test files, silence irrelevant and noisy pylint warnings
# pylint: disable=no-self-use,protected-access,missing-docstring

import unittest

from antismash.common import path, utils
from antismash.modules.nrps_pks.kr_analysis.kr_analysis import run_kr_analysis


class TestKRAnalysis(unittest.TestCase):
    def setUp(self):
        self.query_data = utils.read_fasta(path.get_full_path(__file__, 'data', 'SCO_genes.fasta'))

    def test_full_run(self):
        active, stereo = run_kr_analysis(self.query_data)
        assert active == {'SCO0111_KR1': False,
                          'SCO0126_KR1': True,
                          'SCO0256_KR1': False,
                          'SCO5086_KR1': False,
                          'SCO5097_KR1': False,
                          'SCO6249_KR1': False,
                          'SCO6264_KR1': False,
                          'SCO6273_KR1': True,
                          'SCO6274_KR1': True,
                          'SCO6274_KR2': True,
                          'SCO6275_KR1': True,
                          'SCO6275_KR2': True,
                          'SCO6282_KR1': False,
                          'SCO6827_KR1': True,
                          'SCO6829_KR1': False}
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
