# License: GNU Affero General Public License v3 or later
# A copy of GNU AGPL v3 should have been included in this software package in LICENSE.txt.

# for test files, silence irrelevant and noisy pylint warnings
# pylint: disable=use-implicit-booleaness-not-comparison,protected-access,missing-docstring,invalid-name,too-few-public-methods

import string
import unittest
from unittest.mock import patch

from Bio import SearchIO

from antismash.common import subprocessing
from antismash.common import fasta, path, secmet
from antismash.common.test.helpers import DummyAntismashDomain, DummyRecord, DummyPFAMDomain
from antismash.detection.nrps_pks_domains.modular_domain import TOOL
from antismash.modules import active_site_finder
from antismash.modules.active_site_finder import analysis


def parse_hmmpfam_results(filename):
    with open(path.get_full_path(__file__, 'data', filename), encoding="utf-8") as handle:
        return list(SearchIO.parse(handle, "hmmer2-text"))


def rebuild_domains(filename, domain_type):
    full_path = path.get_full_path(__file__, 'data', filename)
    domain_fasta = fasta.read_fasta(full_path)
    domains = []
    for name, translation in domain_fasta.items():
        domain = DummyAntismashDomain(start=1, end=100, domain_id=domain_type + name,
                                      tool=TOOL)
        domain.domain = domain_type
        domain.translation = translation
        domains.append(domain)
    return domains


class DummyAlignment:
    def __init__(self, domain, positions, scaffold_match=True):
        self.domain = domain
        self.positions = positions
        self.scaffold_match = scaffold_match

    def get_signature(self, positions):
        return "".join([self.positions[position] for position in positions])


class TestAnalyses(unittest.TestCase):
    def setUp(self):
        self.record = secmet.Record()
        # except for Thioesterase, all domains were found in BN001301.1
        # TE domains were found in Y16952
        for filename, domain_type in [("PKS_KS.input", "PKS_KS"), ("AT.input", "PKS_AT"),
                                      ("ACP.input", "ACP"), ("DH.input", "PKS_DH"),
                                      ("KR.input", "PKS_KR"), ("TE.input", "Thioesterase"),
                                      ("ER.input", "PKS_ER")]:
            for domain in rebuild_domains(filename, domain_type):
                self.record.add_antismash_domain(domain)
        # these PFAMs found in BN001301.1 with clusterhmmer, one was excluded
        # to avoid a Biopython SearchIO bug
        domain_fasta = fasta.read_fasta(path.get_full_path(__file__, 'data', "p450.input"))
        for name, translation in domain_fasta.items():
            pfam_domain = DummyPFAMDomain(domain="p450", domain_id="PFAM_p450_" + name)
            pfam_domain.translation = translation
            self.record.add_pfam_domain(pfam_domain)

    def test_KS_N(self):
        with patch.object(subprocessing, "run_hmmpfam2", return_value=parse_hmmpfam_results("KS_N.output")):
            results = analysis.asp_ks(self.record)
        expected = {"PKS_KS0": 'found active site cysteine: True, scaffold matched GSSS: True',
                    "PKS_KS1": 'found active site cysteine: False, scaffold matched GSSS: False',
                    "PKS_KS2": 'found active site cysteine: True, scaffold matched GSSS: False',
                    "PKS_KS3": 'found active site cysteine: True, scaffold matched GSSS: False'}
        assert {dom.domain_id: message for dom, message in results} == expected

    def test_KS_C(self):
        with patch.object(subprocessing, "run_hmmpfam2", return_value=parse_hmmpfam_results("KS_C.output")):
            results = analysis.asp_ks_c(self.record)
        expected = {"PKS_KS0": 'found active site histidines: True',
                    "PKS_KS1": 'found active site histidines: True',
                    "PKS_KS2": 'found active site histidines: True',
                    "PKS_KS3": 'found active site histidines: True'}
        assert {dom.domain_id: message for dom, message in results} == expected

    def test_AT(self):
        with patch.object(subprocessing, "run_hmmpfam2", return_value=parse_hmmpfam_results("AT.output")):
            results = analysis.asp_at(self.record)
        expected = {"PKS_AT0": "found active site serine: True, scaffold matched GHGE: False",
                    "PKS_AT1": "found active site serine: True, scaffold matched GHGE: True",
                    "PKS_AT2": "found active site serine: True, scaffold matched GHGE: True",
                    "PKS_AT3": "found active site serine: True, scaffold matched GHGE: False"}
        assert {dom.domain_id: message for dom, message in results} == expected

    def test_ACP_type(self):
        with patch.object(subprocessing, "run_hmmpfam2", return_value=parse_hmmpfam_results("ACP.output")):
            results = analysis.acp_type(self.record)
        expected = {'ACP0': 'non-beta-branching ACP',
                    'ACP2': 'non-beta-branching ACP',
                    'ACP3': 'non-beta-branching ACP'}
        assert {dom.domain_id: message for dom, message in results} == expected

    def test_ACP(self):
        with patch.object(subprocessing, "run_hmmpfam2", return_value=parse_hmmpfam_results("ACP.output")):
            results = analysis.asp_acp(self.record)
        expected = {'ACP0': 'found active site serine: True',
                    'ACP2': 'found active site serine: True',
                    'ACP3': 'found active site serine: True'}
        assert {dom.domain_id: message for dom, message in results} == expected

    def test_DH(self):
        with patch.object(subprocessing, "run_hmmpfam2", return_value=parse_hmmpfam_results("DH.output")):
            results = analysis.asp_pksi_dh(self.record)
        expected = {'PKS_DH0': 'catalytic triad H,G,P inconclusive',
                    'PKS_DH1': 'catalytic triad H,G,P inconclusive'}
        assert {dom.domain_id: message for dom, message in results} == expected

    def test_KR_stereo(self):
        with patch.object(subprocessing, "run_hmmpfam2", return_value=parse_hmmpfam_results("KR.output")):
            results = analysis.pksi_kr_stereo(self.record)
        expected = {'PKS_KR1': 'KR domain putatively catalyzing D-configuration product formation',
                    'PKS_KR2': 'KR domain putatively catalyzing D-configuration product formation'}
        assert {dom.domain_id: message for dom, message in results} == expected

    def test_KR(self):
        with patch.object(subprocessing, "run_hmmpfam2", return_value=parse_hmmpfam_results("KR.output")):
            results = analysis.asp_pksi_kr(self.record)
        expected = {'PKS_KR1': 'catalytic triad S,Y,N inconclusive',
                    'PKS_KR2': 'catalytic triad S,Y,N found: True'}
        assert {dom.domain_id: message for dom, message in results} == expected

    def test_thioesterase(self):
        with patch.object(subprocessing, "run_hmmpfam2", return_value=parse_hmmpfam_results("TE.output")):
            results = analysis.asp_thioesterase(self.record)
        expected = {'Thioesterase0': 'active site serine present: True'}
        assert {dom.domain_id: message for dom, message in results} == expected

    def test_ER(self):
        with patch.object(subprocessing, "run_hmmpfam2", return_value=parse_hmmpfam_results("ER.output")):
            results = analysis.pksi_er_stereo(self.record)
        expected = {'Thioesterase0': 'active site serine present: True'}
        expected = {'PKS_ER0': 'ER configuration inconclusive',
                    'PKS_ER1': 'ER domain putatively catalyzing 2R-configuration product formation',
                    'PKS_ER2': 'ER configuration inconclusive'}
        assert {dom.domain_id: message for dom, message in results} == expected

    def test_AT_spec(self):
        with patch.object(subprocessing, "run_hmmpfam2", return_value=parse_hmmpfam_results("AT.output")):
            results = analysis.pksi_at_spec(self.record)
        expected = {'PKS_AT0': '(Methyl)Malonyl-CoA specificity inconclusive',
                    'PKS_AT1': 'Neither malonyl-CoA or methylmalonyl-CoA specific',
                    'PKS_AT2': '(Methyl)Malonyl-CoA specificity inconclusive',
                    'PKS_AT3': '(Methyl)Malonyl-CoA specificity inconclusive'}
        assert {dom.domain_id: message for dom, message in results} == expected

    def test_p450(self):
        with patch.object(subprocessing, "run_hmmpfam2", return_value=parse_hmmpfam_results("p450.output")):
            results = analysis.asp_p450_oxy(self.record)
        expected = {'PFAM_p450_0': 'active site cysteine inconclusive',
                    'PFAM_p450_1': 'active site cysteine inconclusive'}
        assert {dom.domain_id: message for dom, message in results} == expected


class TestSynthetic(unittest.TestCase):
    """ Covers those areas that the real inputs above miss """
    def setUp(self):
        self.record = DummyRecord()

    @patch.object(active_site_finder.common.ActiveSiteAnalysis, "get_alignments")
    def test_kr_stereo(self, mocked_get):
        for char in string.ascii_uppercase:
            mocked_get.return_value = [DummyAlignment("match", {102: char})]
            pairings = analysis.pksi_kr_stereo(self.record)
            if char == "D":
                assert pairings == [("match", "KR domain putatively catalyzing D-configuration product formation")]
            else:
                assert pairings == [("match", "KR domain putatively catalyzing L-configuration product formation")]

    @patch.object(active_site_finder.common.ActiveSiteAnalysis, "get_alignments")
    def test_acp_type(self, mocked_get):
        for current in string.ascii_uppercase:
            mocked_get.return_value = [DummyAlignment(current, {37: current})]
            pairings = analysis.acp_type(self.record)
            if current == "W":
                assert pairings == [(current, "beta-branching ACP")]
            else:
                assert pairings == [(current, "non-beta-branching ACP")]

    @patch.object(active_site_finder.common.ActiveSiteAnalysis, "scaffold_matches")
    @patch.object(active_site_finder.common.ActiveSiteAnalysis, "get_alignments")
    def test_pksi_dh(self, mocked_get, mocked_match):
        mocked_match.return_value = True
        for values in ["GGP", "HPP", "HGG", "NNN", "HGP"]:
            mocked_get.return_value = [DummyAlignment(values, {5: values[0], 39: values[1], 44: values[2]})]
            pairings = analysis.asp_pksi_dh(self.record)
            assert pairings == [(values, f"catalytic triad H,G,P found: {values == 'HGP'}")]
        # and an inconclusive
        mocked_match.return_value = False
        # values here don't matter
        mocked_get.return_value = [DummyAlignment("incon", {})]
        pairings = analysis.asp_pksi_dh(self.record)
        assert pairings == [("incon", "catalytic triad H,G,P inconclusive")]

    @patch.object(active_site_finder.common.ActiveSiteAnalysis, "scaffold_matches")
    @patch.object(active_site_finder.common.ActiveSiteAnalysis, "get_alignments")
    def test_thioesterase(self, mocked_get, mocked_match):
        mocked_match.return_value = True
        for char in string.ascii_uppercase:
            mocked_get.return_value = [DummyAlignment(char, {81: char})]
            pairings = analysis.asp_thioesterase(self.record)
            assert pairings == [(char, f"active site serine present: {char == 'S'}")]
        # and an inconclusive
        mocked_match.return_value = False
        # values here don't matter
        mocked_get.return_value = [DummyAlignment("incon", {})]
        pairings = analysis.asp_thioesterase(self.record)
        assert pairings == [("incon", "active site serine inconclusive")]

    @patch.object(active_site_finder.common.ActiveSiteAnalysis, "scaffold_matches")
    @patch.object(active_site_finder.common.ActiveSiteAnalysis, "get_alignments")
    def test_er_stereo(self, mocked_get, mocked_match):
        mocked_match.return_value = True
        for char in string.ascii_uppercase:
            mocked_get.return_value = [DummyAlignment(char, {39: char})]
            pairings = analysis.pksi_er_stereo(self.record)
            if char == "Y":
                assert pairings == [(char, "ER domain putatively catalyzing 2S-configuration product formation")]
            else:
                assert pairings == [(char, "ER domain putatively catalyzing 2R-configuration product formation")]
        # and an inconclusive
        mocked_match.return_value = False
        # values here don't matter
        mocked_get.return_value = [DummyAlignment("incon", {})]
        pairings = analysis.pksi_er_stereo(self.record)
        assert pairings == [("incon", "ER configuration inconclusive")]

    @patch.object(active_site_finder.common.ActiveSiteAnalysis, "scaffold_matches")
    @patch.object(active_site_finder.common.ActiveSiteAnalysis, "get_alignments")
    def test_at_spec(self, mocked_get, mocked_match):
        mocked_match.return_value = True
        for sig in ["HF", "YS", "HS", "YF", "NN"]:
            mocked_get.return_value = [DummyAlignment(sig, {195: sig[0], 197: sig[1]})]
            pairings = analysis.pksi_at_spec(self.record)
            if sig == "HF":
                assert pairings == [(sig, "Malonyl-CoA specific")]
            elif sig == "YS":
                assert pairings == [(sig, "Methylmalonyl-CoA specific")]
            else:
                assert pairings == [(sig, "Neither malonyl-CoA or methylmalonyl-CoA specific")]
        # and an inconclusive
        mocked_match.return_value = False
        # values here don't matter
        mocked_get.return_value = [DummyAlignment("incon", {})]
        pairings = analysis.pksi_at_spec(self.record)
        assert pairings == [("incon", "(Methyl)Malonyl-CoA specificity inconclusive")]

    @patch.object(active_site_finder.common.ActiveSiteAnalysis, "scaffold_matches")
    @patch.object(active_site_finder.common.ActiveSiteAnalysis, "get_alignments")
    def test_p450(self, mocked_get, mocked_match):
        mocked_match.return_value = True
        for char in string.ascii_uppercase:
            mocked_get.return_value = [DummyAlignment(char, {407: char})]
            pairings = analysis.asp_p450_oxy(self.record)
            assert pairings == [(char, f"active site cysteine present: {char == 'C'}")]
        # and an inconclusive
        mocked_match.return_value = False
        # values here don't matter
        mocked_get.return_value = [DummyAlignment("incon", {})]
        pairings = analysis.asp_p450_oxy(self.record)
        assert pairings == [("incon", "active site cysteine inconclusive")]
