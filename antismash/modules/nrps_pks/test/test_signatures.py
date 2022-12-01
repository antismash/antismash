# License: GNU Affero General Public License v3 or later
# A copy of GNU AGPL v3 should have been included in this software package in LICENSE.txt.

# for test files, silence irrelevant and noisy pylint warnings
# pylint: disable=use-implicit-booleaness-not-comparison,protected-access,missing-docstring

import unittest
from unittest.mock import patch

from antismash.common import path, fasta, subprocessing
from antismash.common.test.helpers import DummyAntismashDomain
from antismash.modules.nrps_pks import signatures


class TestAngstromGeneration(unittest.TestCase):
    def test_angstrom(self):
        aligns = fasta.read_fasta(path.get_full_path(__file__, 'data', 'nrpspred_aligns.fasta'))
        domain = DummyAntismashDomain(domain_id="query")
        domain.translation = aligns[domain.domain_id].replace("-", "")
        with patch.object(subprocessing, "run_muscle_single", return_value=aligns):
            sig_10aa, sig_34aa  = signatures.get_a_dom_signatures(domain)
        assert sig_10aa == "DAFYLGMMCK"
        assert sig_34aa == "L--SFDASLFEMYLLTGGDRNMYGPTEATMCATW"


class TestHelpers(unittest.TestCase):
    def test_good_sequence(self):
        check = signatures.verify_good_sequence

        assert check("LSFDASLFEMYLLTGGDRNMYGPTEATMCATW")
        assert not check("!LSFD")
        assert not check("LSFD-LL")
        assert not check("{LS")


    def test_extract(self):
        extract = signatures.extract
        sequence = "MAGIC-HAT"

        assert "MT" == extract(sequence, [0, 8])
        assert "IC" == extract(sequence, [3, 5])
        assert "CH" == extract(sequence, [4, 5])
        assert "C-H" == extract(sequence, [4, 5, 6])
