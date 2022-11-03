# License: GNU Affero General Public License v3 or later
# A copy of GNU AGPL v3 should have been included in this software package in LICENSE.txt.

# for test files, silence irrelevant and noisy pylint warnings
# pylint: disable=use-implicit-booleaness-not-comparison,protected-access,missing-docstring

import unittest

from antismash.modules.nrps_pks.minowa.minowa_at import run_minowa_at
from antismash.modules.nrps_pks.minowa.minowa_cal import run_minowa_cal
from antismash.common import path, fasta
from antismash.config import build_config, destroy_config


class TestMinowaAT(unittest.TestCase):
    query_data = fasta.read_fasta(path.get_full_path(__file__, "data", "SCO.fasta"))

    def setUp(self):
        build_config([])

    def tearDown(self):
        destroy_config()

    def test_full_run(self):
        results = run_minowa_at(self.query_data)
        assert len(results) == len(self.query_data)
        assert set(results) == set(self.query_data)
        results = {key: val.predictions for key, val in results.items()}
        assert results == {'SCO0126_AT1': [('Malonyl-CoA', 81.1),
                                           ('Methoxymalonyl-CoA', 30.9),
                                           ('Methylmalonyl-CoA', 25.6),
                                           ('inactive', 23.2),
                                           ('Propionyl-CoA', 13.8),
                                           ('2-Methylbutyryl-CoA', 12.2),
                                           ('fatty_acid', 7.9),
                                           ('Isobutyryl-CoA', 1.8),
                                           ('CHC-CoA', 1.1),
                                           ('trans-1,2-CPDA', 0.0),
                                           ('Benzoyl-CoA', 0.0),
                                           ('Acetyl-CoA', 0.0),
                                           ('3-Methylbutyryl-CoA', 0.0),
                                           ('Ethylmalonyl-CoA', -3.2)],
                           'SCO0127_AT1': [('Methoxymalonyl-CoA', 29.2),
                                           ('Methylmalonyl-CoA', 26.5),
                                           ('Malonyl-CoA', 22.1),
                                           ('Ethylmalonyl-CoA', 13.7),
                                           ('trans-1,2-CPDA', 0.0),
                                           ('inactive', 0.0),
                                           ('fatty_acid', 0.0),
                                           ('Isobutyryl-CoA', 0.0),
                                           ('CHC-CoA', 0.0),
                                           ('Benzoyl-CoA', 0.0),
                                           ('Acetyl-CoA', 0.0),
                                           ('3-Methylbutyryl-CoA', 0.0),
                                           ('Propionyl-CoA', -0.2),
                                           ('2-Methylbutyryl-CoA', -4.3)],
                           'SCO5892_AT1': [('Malonyl-CoA', 151.7),
                                           ('inactive', 95.9),
                                           ('Methoxymalonyl-CoA', 74.7),
                                           ('Methylmalonyl-CoA', 70.4),
                                           ('Ethylmalonyl-CoA', 43.0),
                                           ('Propionyl-CoA', 35.7),
                                           ('Isobutyryl-CoA', 31.9),
                                           ('CHC-CoA', 27.7),
                                           ('2-Methylbutyryl-CoA', 26.1),
                                           ('Benzoyl-CoA', 25.0),
                                           ('Acetyl-CoA', 13.9),
                                           ('trans-1,2-CPDA', 13.7),
                                           ('3-Methylbutyryl-CoA', 12.5),
                                           ('fatty_acid', 9.7)],
                           'SCO6273_AT1': [('Malonyl-CoA', 171.9),
                                           ('inactive', 73.8),
                                           ('Methoxymalonyl-CoA', 62.1),
                                           ('Methylmalonyl-CoA', 40.8),
                                           ('Propionyl-CoA', 29.3),
                                           ('Acetyl-CoA', 18.6),
                                           ('Isobutyryl-CoA', 15.6),
                                           ('2-Methylbutyryl-CoA', 14.1),
                                           ('Benzoyl-CoA', 9.6),
                                           ('trans-1,2-CPDA', 0.0),
                                           ('fatty_acid', 0.0),
                                           ('Ethylmalonyl-CoA', 0.0),
                                           ('CHC-CoA', 0.0),
                                           ('3-Methylbutyryl-CoA', 0.0)],
                           'SCO6274_AT1': [('Malonyl-CoA', 171.9),
                                           ('inactive', 73.8),
                                           ('Methoxymalonyl-CoA', 62.1),
                                           ('Methylmalonyl-CoA', 40.8),
                                           ('Propionyl-CoA', 29.3),
                                           ('Acetyl-CoA', 18.6),
                                           ('Isobutyryl-CoA', 15.6),
                                           ('2-Methylbutyryl-CoA', 14.1),
                                           ('Benzoyl-CoA', 9.6),
                                           ('trans-1,2-CPDA', 0.0),
                                           ('fatty_acid', 0.0),
                                           ('Ethylmalonyl-CoA', 0.0),
                                           ('CHC-CoA', 0.0),
                                           ('3-Methylbutyryl-CoA', 0.0)],
                           'SCO6274_AT2': [('Malonyl-CoA', 171.9),
                                           ('inactive', 73.8),
                                           ('Methoxymalonyl-CoA', 62.1),
                                           ('Methylmalonyl-CoA', 40.8),
                                           ('Propionyl-CoA', 29.3),
                                           ('Acetyl-CoA', 18.6),
                                           ('Isobutyryl-CoA', 15.6),
                                           ('2-Methylbutyryl-CoA', 14.1),
                                           ('Benzoyl-CoA', 9.6),
                                           ('trans-1,2-CPDA', 0.0),
                                           ('fatty_acid', 0.0),
                                           ('Ethylmalonyl-CoA', 0.0),
                                           ('CHC-CoA', 0.0),
                                           ('3-Methylbutyryl-CoA', 0.0)],
                           'SCO6275_AT1': [('Malonyl-CoA', 209.2),
                                           ('inactive', 103.5),
                                           ('Methoxymalonyl-CoA', 75.4),
                                           ('Methylmalonyl-CoA', 68.4),
                                           ('Isobutyryl-CoA', 37.8),
                                           ('2-Methylbutyryl-CoA', 31.3),
                                           ('Benzoyl-CoA', 30.9),
                                           ('Acetyl-CoA', 30.9),
                                           ('Propionyl-CoA', 29.8),
                                           ('Ethylmalonyl-CoA', 28.1),
                                           ('fatty_acid', 20.5),
                                           ('CHC-CoA', 16.6),
                                           ('3-Methylbutyryl-CoA', 15.4),
                                           ('trans-1,2-CPDA', 15.0)],
                           'SCO6275_AT2': [('Malonyl-CoA', 203.5),
                                           ('inactive', 97.1),
                                           ('Methoxymalonyl-CoA', 72.9),
                                           ('Methylmalonyl-CoA', 61.7),
                                           ('Isobutyryl-CoA', 41.7),
                                           ('Propionyl-CoA', 30.9),
                                           ('Ethylmalonyl-CoA', 16.8),
                                           ('Acetyl-CoA', 16.8),
                                           ('2-Methylbutyryl-CoA', 14.2),
                                           ('Benzoyl-CoA', 13.3),
                                           ('3-Methylbutyryl-CoA', 9.0),
                                           ('fatty_acid', 8.4),
                                           ('CHC-CoA', 3.9),
                                           ('trans-1,2-CPDA', 0.0)],
                           'SCO6275_AT3': [('Malonyl-CoA', 207.6),
                                           ('inactive', 105.9),
                                           ('Methoxymalonyl-CoA', 62.0),
                                           ('Methylmalonyl-CoA', 50.9),
                                           ('Propionyl-CoA', 30.8),
                                           ('Ethylmalonyl-CoA', 17.7),
                                           ('Isobutyryl-CoA', 16.7),
                                           ('2-Methylbutyryl-CoA', 15.7),
                                           ('Acetyl-CoA', 15.4),
                                           ('Benzoyl-CoA', 11.6),
                                           ('CHC-CoA', 9.5),
                                           ('trans-1,2-CPDA', 0.0),
                                           ('fatty_acid', 0.0),
                                           ('3-Methylbutyryl-CoA', 0.0)],
                           'SCO6827_AT1': [('Methylmalonyl-CoA', 165.7),
                                           ('Ethylmalonyl-CoA', 150.9),
                                           ('Methoxymalonyl-CoA', 141.2),
                                           ('2-Methylbutyryl-CoA', 118.3),
                                           ('Malonyl-CoA', 106.6),
                                           ('trans-1,2-CPDA', 94.3),
                                           ('Benzoyl-CoA', 90.8),
                                           ('Isobutyryl-CoA', 90.1),
                                           ('Propionyl-CoA', 89.7),
                                           ('CHC-CoA', 65.8),
                                           ('Acetyl-CoA', 62.2),
                                           ('inactive', 45.4),
                                           ('3-Methylbutyryl-CoA', 43.8),
                                           ('fatty_acid', 23.7)]}


class TestMinowaCAL(unittest.TestCase):
    query_data = {'SCO5892_CAL1': ('TYRELDLHAREVATHLRHAGVGQGPVLLLHPPGLDYLAAFFGC'
                                   'LYAGAVAVPAYPPDNARFGQTVPRLAAIARDCAATHALTTRRV'
                                   'RETVAADGTGRVGTELDGVRWLVTEDLYTGESTTAWEDPGATA'
                                   'RSLAFLQYTSGSTAAPKGVMVEHGNLVRNLRSIHLRLGHDADS'
                                   'GMVSWLPPYHDMGLIGGILTPVYGGFPAHLMAPMTFVQRPLLW'
                                   'LETLSRTGASTSVAPNFGFEQCLRRITPAQRAGLDLSRWRLAL'
                                   'NGAEPIRPDTLDRFAEYFAPAGFDRTALLPCYGLAEATLMVTG'
                                   'VRPADPPVVESFDAAALEAGTARPADPGGVRTTRVVGCGAPVA'
                                   'DVEVAVVDAATGRRVPDGTVAEIRVSGPGVARGYWGRPEAAAE'
                                   'VFGTRIDGEPGNAWLRTGDVGFRHDGQLYVVGRTKDVIIVQGR'
                                   'NIHPQDVEQTAERVGAGL')}

    def setUp(self):
        build_config([])

    def tearDown(self):
        destroy_config()

    def test_full_run(self):
        expected = [("NH2", 118.5), ("fatty_acid", 37.5), ("AHBA", 6.5), ("shikimic_acid", 0.0), ("Acetyl-CoA", -1.6)]
        results = run_minowa_cal(self.query_data)
        assert len(results) == 1
        assert results['SCO5892_CAL1'].predictions == expected
