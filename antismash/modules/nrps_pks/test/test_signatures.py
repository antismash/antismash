# License: GNU Affero General Public License v3 or later
# A copy of GNU AGPL v3 should have been included in this software package in LICENSE.txt.

# for test files, silence irrelevant and noisy pylint warnings
# pylint: disable=use-implicit-booleaness-not-comparison,protected-access,missing-docstring

import unittest
from unittest.mock import Mock, patch

from antismash.common import subprocessing
from antismash.common.test.helpers import DummyAntismashDomain
from antismash.modules.nrps_pks import signatures


class TestAngstromGeneration(unittest.TestCase):
    def test_angstrom(self):
        # bpsC from Y16952, balhimycin
        profile = (
            "KGVmveHrnvvnlvkwlneryflfgeeddllgesdrvLqfssAysFDaSvweifgaLLnGgtLViv"
            "pkefsetrlDpeaLaalieregiTvlnltPsllnllldaaeeatpdfapedlssLrrvlvGGEaLs"
            "pslarrlrerfpdragvrliNaYGPTEtTVcaTi"
        )
        query = (
            "KGVAIPHGAVAGLAGDAG---WQIGPGD-------GVLMHAT-HVFDPSLYAMWVPLVSGAR-VLL"
            "TEP---GVLDAAGVRQAVHR-GATFVHLTAGTFRALAETA--------PECFEGLVEIGTGGDVVP"
            "LQSVENLRRAQ---PGLRVRNTYGPTETTLCATW"
        )
        domain = DummyAntismashDomain(domain_id="query")
        domain.translation = query.replace("-", "")
        hits = [Mock()]
        hsp = Mock()
        hsp.aln = (Mock(), Mock())
        hsp.aln[1].seq = profile
        hsp.aln[0].seq = query
        hsp.hit_start = 0
        hsp.hit_id = signatures.ACTIVE_SITE_PROFILE_NAME
        hits[0].hsps = [hsp]
        with patch.object(subprocessing.hmmpfam, "run_hmmpfam2", return_value=hits) as patched:
            sig_10aa, sig_34aa = signatures.get_a_dom_signatures(domain)
            assert patched.called_once
        assert sig_10aa == "DPYHGGTLCK"
        assert sig_34aa == "LDAAFDPSLYAVHLGTGGDRNTYGPTETTLCATW"


class TestHelpers(unittest.TestCase):
    def test_good_sequence(self):
        check = signatures.verify_good_sequence

        assert check("LSFDASLFEMYLLTGGDRNMYGPTEATMCATW")
        assert not check("!LSFD")
        assert not check("LSFD-LL")
        assert not check("{LS")
