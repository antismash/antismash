# License: GNU Affero General Public License v3 or later
# A copy of GNU AGPL v3 should have been included in this software package in LICENSE.txt.

# for test files, silence irrelevant and noisy pylint warnings
# pylint: disable=no-self-use,protected-access,missing-docstring,invalid-name

import unittest

from Bio import SearchIO
from minimock import mock, restore

from antismash.common import subprocessing  # mocked, pylint: disable=unused-import
from antismash.common import fasta, path, secmet
from antismash.modules.active_site_finder import analysis


def parse_hmmpfam_results(filename):
    with open(path.get_full_path(__file__, 'data', filename)) as handle:
        return list(SearchIO.parse(handle, "hmmer2-text"))


def rebuild_domains(filename, domain_type):
    full_path = path.get_full_path(__file__, 'data', filename)
    domain_fasta = fasta.read_fasta(full_path)
    dummy_location = secmet.feature.FeatureLocation(1, 100)
    domains = []
    for name, translation in domain_fasta.items():
        domain = secmet.feature.AntismashDomain(dummy_location)
        domain.domain = domain_type
        domain.domain_id = domain_type + name
        domain.translation = translation
        domains.append(domain)
    return domains


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
        dummy_location = secmet.feature.FeatureLocation(1, 100)
        domain_fasta = fasta.read_fasta(path.get_full_path(__file__, 'data', "p450.input"))
        for name, translation in domain_fasta.items():
            pfam_domain = secmet.feature.PFAMDomain(dummy_location, description="test")
            pfam_domain.translation = translation
            pfam_domain.domain_id = "PFAM_p450_" + name
            pfam_domain.domain = "p450"
            self.record.add_pfam_domain(pfam_domain)

    def tearDown(self):
        restore()

    def test_KS_N(self):
        mock("subprocessing.run_hmmpfam2", returns=parse_hmmpfam_results("KS_N.output"))
        results = analysis.asp_ks(self.record)
        expected = {"PKS_KS0": 'found active site cysteine: False, scaffold matched GSSS: True',
                    "PKS_KS1": 'found active site cysteine: False, scaffold matched GSSS: False',
                    "PKS_KS2": 'found active site cysteine: False, scaffold matched GSSS: False',
                    "PKS_KS3": 'found active site cysteine: False, scaffold matched GSSS: False'}
        assert {dom.domain_id: message for dom, message in results} == expected

    def test_KS_C(self):
        mock("subprocessing.run_hmmpfam2", returns=parse_hmmpfam_results("KS_C.output"))
        results = analysis.asp_ks_c(self.record)
        expected = {"PKS_KS0": 'found active site histindines: False',
                    "PKS_KS1": 'found active site histindines: False',
                    "PKS_KS2": 'found active site histindines: False',
                    "PKS_KS3": 'found active site histindines: False'}
        assert {dom.domain_id: message for dom, message in results} == expected

    def test_AT(self):
        mock("subprocessing.run_hmmpfam2", returns=parse_hmmpfam_results("AT.output"))
        results = analysis.asp_at(self.record)
        expected = {"PKS_AT0": "found active site serine: True, scaffold matched GHGE: False",
                    "PKS_AT1": "found active site serine: False, scaffold matched GHGE: True",
                    "PKS_AT2": "found active site serine: False, scaffold matched GHGE: True",
                    "PKS_AT3": "found active site serine: False, scaffold matched GHGE: False"}
        assert {dom.domain_id: message for dom, message in results} == expected

    def test_ACP_type(self):
        mock("subprocessing.run_hmmpfam2", returns=parse_hmmpfam_results("ACP.output"))
        results = analysis.acp_type(self.record)
        expected = {'ACP0': 'non-beta-branching ACP',
                    'ACP2': 'non-beta-branching ACP',
                    'ACP3': 'non-beta-branching ACP'}
        assert {dom.domain_id: message for dom, message in results} == expected

    def test_ACP(self):
        mock("subprocessing.run_hmmpfam2", returns=parse_hmmpfam_results("ACP.output"))
        results = analysis.asp_acp(self.record)
        expected = {'ACP0': 'found active site serine: True',
                    'ACP2': 'found active site serine: True',
                    'ACP3': 'found active site serine: True'}
        assert {dom.domain_id: message for dom, message in results} == expected

    def test_DH(self):
        mock("subprocessing.run_hmmpfam2", returns=parse_hmmpfam_results("DH.output"))
        results = analysis.asp_pksi_dh(self.record)
        expected = {'PKS_DH0': 'catalytic triade H,G,P inconclusive',
                    'PKS_DH1': 'catalytic triade H,G,P inconclusive'}
        assert {dom.domain_id: message for dom, message in results} == expected

    def test_KR_bug(self):
        # currently breaks in BioPython
        try:
            parse_hmmpfam_results("KR.output")
        except ValueError as err:
            if "Sequence lengths do not match" not in str(err):
                raise
            self.skipTest("Biopython bug preventing real output check")
        # and if it works, update the other KR tests to use KR.output instead of KR.output.altered
        self.fail("Unexpected success parsing hmmpfam2 output, biopython bug probably fixed. Update this test.")

    def test_KR_stereo(self):
        # use an altered output to avoid the Biopython bug
        mock("subprocessing.run_hmmpfam2", returns=parse_hmmpfam_results("KR.output.altered"))
        results = analysis.pksi_kr_stereo(self.record)
        expected = {'PKS_KR1': 'KR domain putatively catalyzing D-configuration product formation',
                    'PKS_KR2': 'KR domain putatively catalyzing D-configuration product formation'}
        assert {dom.domain_id: message for dom, message in results} == expected

    def test_KR(self):
        # use an altered output to avoid the Biopython bug
        mock("subprocessing.run_hmmpfam2", returns=parse_hmmpfam_results("KR.output.altered"))
        results = analysis.asp_pksi_kr(self.record)
        expected = {'PKS_KR1': 'catalytic triade S,Y,N inconclusive',
                    'PKS_KR2': 'catalytic triade S,Y,N found: True'}
        assert {dom.domain_id: message for dom, message in results} == expected

    def test_thioesterase(self):
        mock("subprocessing.run_hmmpfam2", returns=parse_hmmpfam_results("TE.output"))
        results = analysis.asp_thioesterase(self.record)
        expected = {'Thioesterase0': 'active site serine present: True'}
        assert {dom.domain_id: message for dom, message in results} == expected

    def test_ER(self):
        mock("subprocessing.run_hmmpfam2", returns=parse_hmmpfam_results("ER.output"))
        results = analysis.pksi_er_stereo(self.record)
        expected = {'Thioesterase0': 'active site serine present: True'}
        expected = {'PKS_ER0': 'ER configuration inconclusive',
                    'PKS_ER1': 'ER domain putatively catalyzing 2R-configuration product formation',
                    'PKS_ER2': 'ER configuration inconclusive'}
        assert {dom.domain_id: message for dom, message in results} == expected

    def test_AT_spec(self):
        mock("subprocessing.run_hmmpfam2", returns=parse_hmmpfam_results("AT.output"))
        results = analysis.pksi_at_spec(self.record)
        expected = {'PKS_AT0': '(Methyl)Malonyl-CoA specificity inconclusive',
                    'PKS_AT1': 'Neither malonyl-CoA or methylmalonyl-CoA specific',
                    'PKS_AT2': '(Methyl)Malonyl-CoA specificity inconclusive',
                    'PKS_AT3': '(Methyl)Malonyl-CoA specificity inconclusive'}
        assert {dom.domain_id: message for dom, message in results} == expected

    def test_p450(self):
        mock("subprocessing.run_hmmpfam2", returns=parse_hmmpfam_results("p450.output"))
        results = analysis.asp_p450_oxy(self.record)
        expected = {'PFAM_p450_0': 'active site cystein inconclusive',
                    'PFAM_p450_1': 'active site cystein inconclusive'}
        assert {dom.domain_id: message for dom, message in results} == expected
