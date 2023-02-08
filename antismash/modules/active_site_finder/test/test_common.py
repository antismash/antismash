# License: GNU Affero General Public License v3 or later
# A copy of GNU AGPL v3 should have been included in this software package in LICENSE.txt.

# for test files, silence irrelevant and noisy pylint warnings
# pylint: disable=use-implicit-booleaness-not-comparison,protected-access,missing-docstring,consider-using-with

import unittest
from unittest.mock import patch

from Bio import SearchIO

from antismash.common import subprocessing
from antismash.common import fasta, path
from antismash.common.secmet.features import FeatureLocation
from antismash.common.test.helpers import DummyAntismashDomain, DummyPFAMDomain
from antismash.modules.active_site_finder.common import ActiveSiteAnalysis, Alignment, get_signature


class TestCommon(unittest.TestCase):
    def test_get_signature(self):
        query = ("TYLVTGGAGGIGGQLALWLAD-QGARHLLLTGRS-A-L--PEQdavvsethpqaTAVAVLRQLRERGVNVTYKAVDVADAH"
                 "AMQATLESRRRA-GM--PPVRGVFHAAGVIDYTLLSDMSGAEMDRVLAAKVSGAWNLHRLLR-EES----VEAFVLFSSGS"
                 "ALLSSPMLGGYAAGNAFLDALAHHRHAQGL--SGTVVNWGFWD--")
        aln = ("tYLitGGlGGLGlslArWLaerrGARrLvLlSRslglpllpsp...........eaqellaeLealGarVrvvacDVtdraav"
               "rrllaeiraldtlespPirGViHaAgVLrDallenmtaedfrrVlaPKVdGawnLHeatreddppegsLDFFvlFSSiagllG"
               "npGQanYAAANaFLDAlAryRRarGLRGpAlsinWGaWadv")
        positions = [102]

        assert get_signature(query, aln, positions) == "Y"

        # bad coords
        assert get_signature(query, aln, [len(aln)]) == ""


class TestAlignment(unittest.TestCase):
    def setUp(self):
        self.domain = DummyPFAMDomain(domain="p450")
        self.alignment = Alignment(self.domain, "WLAD-QGAR", "WLaer.rGA", 10, 19)

    def test_get_signature(self):
        assert self.alignment.get_signature([11, 13, 15]) == "WA-"

    def test_bad_positions(self):
        with self.assertRaisesRegex(ValueError, "Hit start position is negative"):
            Alignment(self.domain, "A", "B", -1, 10)
        with self.assertRaisesRegex(ValueError, "Hit end position is not greater than the start"):
            Alignment(self.domain, "A", "B", 10, 1)


class TestAnalysisCore(unittest.TestCase):
    def generate_domains(self):
        inputs = fasta.read_fasta(path.get_full_path(__file__, 'data', 'PKS_KS.input'))
        domains = []
        last_end = 0
        for translation in inputs.values():
            location = FeatureLocation(last_end + 10, last_end + len(translation)*3 + 16)
            domain = DummyAntismashDomain(location=location)
            domain.translation = translation
            domains.append(domain)
            domain.domain = "PKS_KS"

        location = FeatureLocation(last_end + 10, last_end + len(domains[-1].translation)*3 + 16)
        domains.append(DummyAntismashDomain(location=location))
        domains[-1].domain = "PKS_KR"
        return domains

    def test_bad_args(self):
        with self.assertRaisesRegex(ValueError, "No database file located for"):
            ActiveSiteAnalysis("test", tuple(), "bad_db", [5, 6], ["C", "S"])
        with self.assertRaisesRegex(ValueError, "Number of expected values must match number of positions"):
            ActiveSiteAnalysis("test", tuple(), "PKSI-KR.hmm2", [5, 6], ["C"])
        with self.assertRaisesRegex(ValueError, "Number of expected values must match number of positions"):
            ActiveSiteAnalysis("test", tuple(), "PKSI-KR.hmm2", [5], ["C", "S"])
        with self.assertRaisesRegex(ValueError, "Number of emissions must match number of positions"):
            ActiveSiteAnalysis("test", tuple(), "PKSI-KR.hmm2", [5, 6], ["C", "S"], emissions=[1.0, .2, .3])
        with self.assertRaisesRegex(ValueError, "Number of emissions must match number of positions"):
            ActiveSiteAnalysis("test", tuple(), "PKSI-KR.hmm2", [5, 6], ["C", "S"], emissions=[1.0])

    def test_no_candidates_doesnt_break(self):
        analysis = ActiveSiteAnalysis("test", tuple(), "PKSI-KR.hmm2", [5, 6], ["C", "S"])
        with patch.object(subprocessing, "run_hmmpfam2", side_effect=RuntimeError()) as mocked_run:
            assert analysis.get_alignments() == []
            mocked_run.assert_not_called()

    def test_domains_of_interest(self):
        domain = DummyPFAMDomain(domain="p450")
        analysis = ActiveSiteAnalysis("not-p450", (domain,), "PKSI-KR.hmm2", [5, 6], ["C", "S"])
        assert analysis.domains_of_interest == []
        analysis = ActiveSiteAnalysis("p450", (domain,), "PKSI-KR.hmm2", [5, 6], ["C", "S"])
        assert analysis.domains_of_interest == [domain]

    def test_bad_candidates(self):
        with self.assertRaisesRegex(TypeError, "Candidates must be Domains, not"):
            ActiveSiteAnalysis("not-p450", ("not a Domain",), "PKSI-KR.hmm2", [5, 6], ["C", "S"])

    def test_alignment_generation(self):
        pregenerated = list(SearchIO.parse(open(path.get_full_path(__file__, 'data', 'KS_N.output'),
                                                encoding="utf-8"),
                                           "hmmer2-text"))
        domains = self.generate_domains()
        analysis = ActiveSiteAnalysis("PKS_KS", domains, "PKSI-KS_N.hmm2",
                                      [176, 186, 187, 188], ['G', 'S', 'S', 'S'])
        with patch.object(subprocessing, "run_hmmpfam2", return_value=pregenerated):
            alignments = analysis.get_alignments()
        assert {"PKS_KS"} == {domain.domain for domain in analysis.domains_of_interest}
        assert len(alignments) == 4
        assert [align.domain for align in alignments[:4]] == domains[:4]

    def test_alignment_generation_no_hits(self):
        no_hits = list(SearchIO.parse(open(path.get_full_path(__file__, 'data', 'no_hits.output'),
                                           encoding="utf-8"),
                                      "hmmer2-text"))
        for result in no_hits:
            assert not result.hsps, "hits shouldn't exist"
        domains = self.generate_domains()
        analysis = ActiveSiteAnalysis("PKS_KS", domains, "PKSI-KS_N.hmm2",
                                      [176, 186, 187, 188], ['G', 'S', 'S', 'S'])
        with patch.object(subprocessing, "run_hmmpfam2", return_value=no_hits):
            assert analysis.get_alignments() == []


class TestScaffoldMatching(unittest.TestCase):
    def setUp(self):
        domain = DummyPFAMDomain(domain="p450")
        self.alignment = Alignment(domain, "WLAD-QGAR", "WLae.rGAR", 10, 19)

    def create_analysis(self, positions, expected):
        return ActiveSiteAnalysis("test", tuple(), "PKSI-KR.hmm2", positions, expected)

    def test_mismatch_over_gap(self):
        # WDG vs WEG
        analysis = self.create_analysis([1, 4, 7], ["W", "e", "G"])
        assert not analysis.scaffold_matches(self.alignment)

    def test_mismatch_on_gap(self):
        # W-G vs WrG
        analysis = self.create_analysis([1, 5, 7], ["W", "r", "G"])
        assert not analysis.scaffold_matches(self.alignment)

    def test_match_no_gaps(self):
        # WL
        analysis = self.create_analysis([1, 2], ["W", "L"])
        assert analysis.scaffold_matches(self.alignment)

    def test_match_no_gaps_case(self):
        # WLA vs WLa
        analysis = self.create_analysis([1, 2, 3], ["W", "L", "a"])
        assert not analysis.scaffold_matches(self.alignment)

    def test_match_over_gap(self):
        # WLG
        analysis = self.create_analysis([1, 2, 6], ["W", "L", "G"])
        assert analysis.scaffold_matches(self.alignment)
