# License: GNU Affero General Public License v3 or later
# A copy of GNU AGPL v3 should have been included in this software package in LICENSE.txt.

# for test files, silence irrelevant and noisy pylint warnings
# pylint: disable=use-implicit-booleaness-not-comparison,protected-access,missing-docstring

import os
import unittest

from antismash.common import path
from antismash.common.test.helpers import DummyCDS
from antismash.config import build_config, destroy_config
from antismash.modules.nrps_pks import orderfinder


class TestCTerminalExtract(unittest.TestCase):
    def setUp(self):
        # 110 long, since only the last 100 will be used
        self.seqs = {"STAUR_3982": ("FLEFTRQRGFISEEFGREHDSELMKTYLPTLRKDLVLLESYSYAEEAPLDMPLTV"
                                    "FASTRDRIIPSTQLESWGELTREKPSIHLFEGDHFFARDAGGPLLALIREKLGLG"),
                     "STAUR_3984": ("LARVLRMEASRIDRLRALGELGLDSLMSLELRNRLEASLGLKLSVTLLFTYPNLA"
                                    "GLAEYLHGELLPAAAREQPAAQSQTHAAPSQIAEQVEQLSKDELLAFFDKSFGIA"),
                     "STAUR_3983": ("TNMGLDSLMSLELRNRLEATLGLKLSATLLFTYPNLAALADHLLGKLSSVDEAPA"
                                    "KTAPTAAAPPPPPTLKPQAALPAELDQLGKDELLSLFDESLTESLKRTRMTRTSR"),
                     "STAUR_3985": ("PSKIDRLRALGELGLDSLMSLELRNRLEAALGMKLSATLLFTYPNLASLAQHVVG"
                                    "RMEFPSEATVAPITASPGAVEGQAERLAEVEQMSDDEAEQLLLASLESLSTELLK"),
                     "STAUR_3972": ("SEAALRGSAAGVAYTASKHALIGFTKNTAFMYGAKGVRVNIVAPGPVRTSISGAS"
                                    "RSDHGWSRIAPVMNVLAVPVAESATLAGHILWLMSDEAENINGAVLPSDGGWSTF")
                     }

        self.data_dir = path.get_full_path(os.path.dirname(__file__), "..", "data", "terminals")

        self.features = [DummyCDS(1, 200, translation=seq, locus_tag=seq_id) for seq_id, seq in self.seqs.items()]
        self.features_by_id = {feature.locus_tag: feature for feature in self.features}
        build_config([])

    def tearDown(self):
        destroy_config()

    def test_c_terminals_no_end(self):
        residues = orderfinder.extract_cterminus(self.data_dir, self.features, "")
        assert residues == {'STAUR_3972': 'ES', 'STAUR_3982': 'GK',
                            'STAUR_3983': 'DS', 'STAUR_3984': 'DS',
                            'STAUR_3985': 'DS'}

    def test_c_terminals_with_end(self):
        residues = orderfinder.extract_cterminus(self.data_dir, self.features, self.features_by_id["STAUR_3982"])
        assert residues == {'STAUR_3972': 'ES', 'STAUR_3983': 'DS',
                            'STAUR_3984': 'DS', 'STAUR_3985': 'DS'}


class TestNTerminalExtract(unittest.TestCase):
    def setUp(self):
        # 60 long, only the first 50 should be used
        self.seqs = {"STAUR_3972": "MQPLEGRFAGRTVVVTGAGAGIGHATASRLMREGARVVASDIAQDRLAALEAESPRGALV",
                     "STAUR_3982": "MSQPENEYLSRLRNAVVALREMQQEIDALNHARTEPIAIVGMGCRFPGGASTPEAFWKLL",
                     "STAUR_3985": "MRQAGSPSSPEALQSLVISLVAARTALPVRSIDVREPLSRHGLDSAGAMGLLAELSADLG",
                     "STAUR_3983": "MSVSEADYIARLRKAAITLKEMEGKLGALERARTEPIAIIGMGCRLPGGASTPEAFWKLL",
                     "STAUR_3984": "MNDASSMSTVKRALLAVQEMKARLDAVTRAQTEPIAIIGLGCRLPGGASTPEAFWKLIES"}

        self.data_dir = path.get_full_path(os.path.dirname(__file__), "..", "data", "terminals")

        self.features = [DummyCDS(1, 200, translation=seq, locus_tag=seq_id) for seq_id, seq in self.seqs.items()]
        self.features_by_id = {feature.locus_tag: feature for feature in self.features}
        build_config([])

    def tearDown(self):
        destroy_config()

    def test_n_terminals_no_start(self):
        residues = orderfinder.extract_nterminus(self.data_dir, self.features, None)
        assert residues == {'STAUR_3972': 'L-', 'STAUR_3982': 'ER',
                            'STAUR_3983': 'DK', 'STAUR_3984': 'SQ',
                            'STAUR_3985': 'SV'}

    def test_n_terminals_with_start(self):
        residues = orderfinder.extract_nterminus(self.data_dir, self.features, self.features_by_id["STAUR_3982"])
        assert residues == {'STAUR_3972': 'L-', 'STAUR_3983': 'DK',
                            'STAUR_3984': 'SQ', 'STAUR_3985': 'SV'}
