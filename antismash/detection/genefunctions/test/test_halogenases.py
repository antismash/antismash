# License: GNU Affero General Public License v3 or later
# A copy of GNU AGPL v3 should have been included in this software package in LICENSE.txt.

# for test files, silence irrelevant and noisy pylint warnings
# pylint: disable=use-implicit-booleaness-not-comparison,protected-access,missing-docstring

import json
import random
import unittest
from unittest.mock import patch

from antismash.common import utils
from antismash.common import secmet
from antismash.common import subprocessing
# fake classes for testing
from antismash.common.secmet.test.helpers import (
    DummyRecord,
    DummyFeature,
    DummyCDS)
from antismash.common.test.helpers import FakeHSPHit, FakeHit
# functions and classes to test
from antismash.detection.genefunctions import halogenases
from antismash.detection.genefunctions.halogenases import (
    HalogenaseResult,
    HalogenaseHmmResult,
    Match,
    HalogenaseCategories,
    FlavinDependents,
    TryptophanSubstrate,
    run_halogenase_phmms,
    search_signature_residues,
    check_for_fdh,
    check_for_halogenases,
    specific_analysis
)

test_fasta = {
        "ktzQ": """MDDNRIRSILVLGGGTAGWMSACYLSKALGPGVEVTVLEAPSISRIRVGEATIPNLHKVF
                FDFLGIAEDEWMRECNASYKAAVRFVNWRTPGDGQATPRRRPDGRPDHFDHLFGQLPEHE
                NLPLSQYWAHRRLNGLTDEPFDRSCYVQPELLDRKLSPRLMDGTKLASYAWHFDADLVAD
                FLCRFAVQKLNVTHVQDVFTHADLDQRGHITAVNTESGRTLAADLFIDCSGFRSVLMGKV
                MQEPFLDMSKHLLNDRAVALMLPHDDEKVGIEPYTSSLAMRSGWSWKIPLLGRFGSGYVY
                SSQFTSQDEAAEELCRMWDVDPAEQTFNNVRFRVGRSRRAWVRNCVAIGVSAMFVEPLES
                TGLYFSYASLYQLVKHFPDKRFRPILADRFNREVATMYDDTRDFLQAHFSLSPRDDSEFW
                RACKELPFADGFAEKVEMYRAGLPVELPVTIDDGHYYGNFEAEFRNFWTNSNYYCIFAGL
                GFLPEHPLPVLEFRPEAVDRAEPVFAAVRRRTEELVATAPTMQAYLRRLHQGT""",
        "ktzR": """MTAAYLKTAFGDRLSITVVESSRIGTIGVGEATFSDIQHFFQFL
                NLREQDWMPACNATYKLGIRFENWRHVGHHFYQPFEQIRPVYGFPLTDWWLHDAPTDR
                FDTDCFVMPNLCEAGRSPRHLDGTLADEDFVEEGDELANRTMSEHQGKSQFPYAYHFE
                AALLAKFLTGYAVDRGVEHVVDDVLDVRLDQRGWIEHVVTAEHGEIHGDLFVDCTGFR
                GLLLNKALGVPFVSYQDTLPNDSAVALQVPLDMQRRGIVPNTTATAREAGWIWTIPLF
                GRVGTGYVYAKDYLSPEEAERTLREFVGPAAADVEANHIRMRIGRSQESWRNNCVAIG
                LSSGFVEPLESTGIFFIHHAIEQLVKHFPAADWNPKSRDMYNSAVAHVMDGIREFLVI
                HYRGAARADNQYWRDTKTRPLPDGLAERIECWQTQLPDTETIYPYYHGLPPYSYMCIL
                MGGGAIRTPASAALALTDQGAAQKEFAAVRDRAAQLRDTLPSHYEYLARMRGLDV""",
        "mibH": """MLKNVVVVGGGTAGWMTASYLTAAFGDRIGVTLVESKRVGSIGVGEATFSTVRHFFEYLG
                LEEKEWMPACNATYKLAIRFENWREPGHHFYHPFERQRVVDGFPLTDWWLREPRSDRFDK
                DCFLVGTLCDDLKSPRQLNGELFEGGLGGRSAYRTTLAEQTTQFPYAYHFDATLVANYLR
                DYAVARGVKHVLDDVQDVALDDRGWISHVVTGESGNLTGDLFIDCTGFRSLLLGKALAEP
                FQSYQDSLPNDSAVALRVPQDMENRGLRPCTTATAQEAGWIWTIPLFDRIGTGYVYAGDY
                ISPEEAERTLRAFVGPAAEHADANHIKMRIGRSNRHWVNNCVAVGLSSGFVEPLESTGIF
                FIQHAIEQLVKHFPDERWDDGLRTAYNKLVNNVMDGVREFLVVHYYAAKRQDNQYWKDAK
                TRPLPDGLAERLERWQTRLPDNESVFPHYHGFESYSYVCMLLGLGGLDLKSSPALGLMDA
                APARHEFKLVGEQAAELARTLPTQYEYFAQLHRAR"""}


class TestHalogenasesAnalysis(unittest.TestCase):
    def setUp(self):
        self.test_trp_5_match = Match("trp_5", 5, 0, "")
        self.test_trp_6_7_match = Match("trp_6_7", 6, 1, "")
        self.trp_6_7_hmm_result = HalogenaseHmmResult(
            hit_id='trp_6_7',
            bitscore=1000,
            query_id='trp_6_7',
            enzyme_type='Flavin-dependent',
            profile=halogenases.pHMM_SIGNATURES[1].hmm_file,
        )
        self.trp_5_hmm_result = HalogenaseHmmResult(
            hit_id='trp_5',
            bitscore=1000,
            query_id='trp_5',
            enzyme_type='Flavin-dependent',
            profile=halogenases.pHMM_SIGNATURES[0].hmm_file,
        )
        matches = [self.test_trp_5_match, self.test_trp_6_7_match]
        self.enzyme_with_more_matches = HalogenaseResult("ktzR", potential_matches=matches)

        self.empty_enzyme = HalogenaseResult("ktzQ")

    def test_get_best_match(self):
        assert not self.empty_enzyme.get_best_match()

        positive_test_best_match = self.enzyme_with_more_matches.get_best_match()
        assert len(positive_test_best_match) == 1 and isinstance(positive_test_best_match[0], Match)

        self.enzyme_with_more_matches.add_potential_matches(self.test_trp_6_7_match)
        assert len(self.enzyme_with_more_matches.potential_matches) == 3

        multiple_matches = self.enzyme_with_more_matches.get_best_match()
        assert len(multiple_matches) == 2 and isinstance(positive_test_best_match[0], Match)

        one_potential_match = HalogenaseResult("test_enzyme",
                                               potential_matches=[self.test_trp_5_match])
        multiple_matches = one_potential_match.get_best_match()
        assert len(multiple_matches) == 1 and isinstance(positive_test_best_match[0], Match)

    def test_conversion_methods(self):
        converted_to_json = self.enzyme_with_more_matches.to_json()
        converted_to_json = json.loads(json.dumps(converted_to_json))
        assert isinstance(converted_to_json, dict)

        converted_from_json = HalogenaseResult.from_json(converted_to_json)
        assert isinstance(converted_from_json, HalogenaseResult)
        assert converted_to_json == converted_from_json.to_json()

    @patch.object(subprocessing, "run_hmmsearch",
                  return_value=[FakeHit("start", "end", 1000, "foo")])
    def test_run_halogenase_phmms(self, run_hmmsearch):
        for value in run_hmmsearch.return_value:
            value.hsps = [FakeHSPHit("foo", "foo", bitscore=250)]

        negative_test_halogenase_hmms_by_id = run_halogenase_phmms("")
        assert not negative_test_halogenase_hmms_by_id

        for value in run_hmmsearch.return_value:
            value.hsps = [FakeHSPHit("foo", "foo", bitscore=1000)]

        positive_test_halogenase_hmms_by_id = run_halogenase_phmms("")
        for hit in positive_test_halogenase_hmms_by_id["foo"]:
            assert isinstance(hit, HalogenaseHmmResult)

    def test_search_signature_residues(self):
        positions = [random.randrange(512, 600) for number in range(10)]

        with patch.object(subprocessing.hmmpfam, "run_hmmpfam2") as run_hmmpfam2:
            run_hmmpfam2.return_value = []
            signature_residues = search_signature_residues(test_fasta["ktzR"],
                                                           positions, self.trp_6_7_hmm_result)
            assert signature_residues is None

            run_hmmpfam2.return_value = [FakeHit("start", "end", 1000, "foo")]
            # checking if hit_id != query_id it breaks
            for result in run_hmmpfam2.return_value:
                result.hsps = [FakeHSPHit("foo", hit_id="foo")]
                for hit in result.hsps:
                    hit_query = DummyFeature()
                    hit_profile = DummyFeature()
                    hit_query.seq = "ghjvghjkbln"
                    hit_profile.seq = "xfdhcgkbjlnkml"
                    hit.aln = [hit_profile, hit_query]

            signature_residues = search_signature_residues(list(test_fasta.values())[0],
                                                           positions, self.trp_6_7_hmm_result)
            assert signature_residues is None

            # checking if hit_id == query_id it runs til the reference position extraction
            for result in run_hmmpfam2.return_value:
                result.hsps = [FakeHSPHit("trp_6_7", hit_id="trp_6_7")]
                for hit in result.hsps:
                    hit_query = DummyFeature()
                    hit_profile = DummyFeature()
                    hit_query.seq = test_fasta["ktzR"]
                    hit_profile.seq = "xfdhcgkbjlnkml"
                    hit.aln = [hit_profile, hit_query]

            with patch.object(utils, "extract_by_reference_positions", return_value="dummy"):
                signature_residues = search_signature_residues(test_fasta["ktzR"],
                                                               positions,
                                                               self.trp_6_7_hmm_result)
                assert signature_residues == "dummy"

    def test_false_check_for_match(self):
        false_test = HalogenaseCategories("fake_name", "visil")
        checked_match = false_test.check_for_match(self.empty_enzyme, self.trp_6_7_hmm_result, 5, [300])
        assert checked_match is False

    def test_check_for_fdh(self):
        with patch.object(halogenases, "get_residues",
                          return_value={"trp_5": FlavinDependents.TRP_5_SIGNATURE_RESIDUES}):
            check_for_fdh(DummyCDS(), self.empty_enzyme, self.trp_5_hmm_result)
            assert self.empty_enzyme.potential_matches[0].profile == "trp_5"
            assert self.empty_enzyme.potential_matches[0].position == 5

        with patch.object(halogenases, "get_residues",
                          return_value={"trp_5": TryptophanSubstrate.TRP_5_SIGNATURE_RESIDUES}):
            low_quality_hit = HalogenaseHmmResult("trp_5", 400, "trp_5", "foo", "trp_5_v2")
            check_for_fdh(DummyCDS(), self.empty_enzyme, low_quality_hit)
            assert self.empty_enzyme.potential_matches[1].profile == "trp_5"
            assert self.empty_enzyme.potential_matches[1].confidence == 0.5

        with patch.object(halogenases, "get_residues",
                          return_value={"trp_6_7": FlavinDependents.TRP_6_SIGNATURE_RESIDUES}):
            check_for_fdh(DummyCDS(), self.empty_enzyme, self.trp_6_7_hmm_result)
            assert self.empty_enzyme.potential_matches[2].profile == "trp_6_7"
            assert self.empty_enzyme.potential_matches[2].position == 6

        with patch.object(halogenases, "get_residues", return_value={"trp_6_7": "VISIL"}):
            check_for_fdh(DummyCDS(), self.empty_enzyme, self.trp_6_7_hmm_result)
            assert self.empty_enzyme.potential_matches[3].profile == "trp_6_7"
            assert self.empty_enzyme.potential_matches[3].position == 7

        with patch.object(halogenases, "get_residues",
                          return_value={"trp_6_7": None}):
            assert check_for_fdh(DummyCDS(), self.empty_enzyme,
                                 self.trp_6_7_hmm_result) is None

        with patch.object(halogenases, "get_residues",
                          return_value={"trp_5": TryptophanSubstrate.TRP_5_SIGNATURE_RESIDUES}):
            substrate = TryptophanSubstrate("trp_5", TryptophanSubstrate.TRP_5_SIGNATURE_RESIDUES)
            with patch.object(halogenases, "TryptophanSubstrate",
                              return_value=substrate):
                low_quality_hit = HalogenaseHmmResult("wrong_name", 400, "trp_5", "foo", "trp_5_v2")
                assert not check_for_fdh(DummyCDS(), self.empty_enzyme, low_quality_hit)

    def test_check_for_halogenases(self):
        negative_checked_halogenases = check_for_halogenases(DummyCDS(), [])
        assert negative_checked_halogenases is None

        mock_cds_feature = DummyCDS(locus_tag="foo", translation="")
        positive_checked_halogenase = check_for_halogenases(mock_cds_feature, [self.trp_5_hmm_result])
        assert isinstance(positive_checked_halogenase, HalogenaseResult)
        assert positive_checked_halogenase.cds_name == "foo"

    def test_get_signatures(self):
        assert TryptophanSubstrate.get_signatures() == [FlavinDependents.TRP_5_SIGNATURE,
                                                        FlavinDependents.TRP_6_SIGNATURE]

    def test_get_residues(self):
        mock_trp_6 = DummyCDS(locus_tag="ktzR",
                              translation=test_fasta["ktzR"])
        tryptophan_accepting = TryptophanSubstrate("ktzR", "")
        trp_6_7_residues = tryptophan_accepting.get_residues(mock_trp_6,
                                                             self.trp_6_7_hmm_result)
        assert (isinstance(trp_6_7_residues, dict) and
                trp_6_7_residues['trp_6_7'] == "TEGCAGFDAYHDRFGNADYGLSIIAKIL")

        mock_trp_5 = DummyCDS(locus_tag="mibH",
                              translation=test_fasta["mibH"])
        tryptophan_accepting = TryptophanSubstrate("mibH", "")
        trp_5_residues = tryptophan_accepting.get_residues(mock_trp_5,
                                                           self.trp_5_hmm_result)
        assert (isinstance(trp_6_7_residues, dict) and
                trp_5_residues['trp_5'] == "VSILIREPGLPRGVPRAVLPGEA")


@patch.object(secmet.Record, "get_cds_by_name", return_value="ktzR")
class TestSpecificAnalysis(unittest.TestCase):
    def setUp(self):
        self.record = DummyRecord()
        self.cds = DummyCDS(locus_tag="ktzR", translation=test_fasta["ktzR"])
        self.test_trp_5_match = Match("trp_5", 5, 0, "")
        self.test_trp_6_7_match = Match("trp_6_7", 6, 1, "")
        self.empty_match = HalogenaseResult("test_enzyme")

    # one best match
    @patch.object(subprocessing, "run_hmmsearch",
                  return_value=[FakeHit("start", "end", 900, "ktzR")])
    @patch.object(secmet.Record, "get_cds_features_within_regions",
                  return_value=[DummyCDS(locus_tag="ktzR", translation=test_fasta["ktzR"])])
    def test_one_best_match(self, _patched_get_cds, run_hmmsearch, _patched_by_name):
        for value in run_hmmsearch.return_value:
            value.hsps = [FakeHSPHit("trp_6_7", "trp_6_7", bitscore=1500)]
        with patch.object(halogenases, "check_for_halogenases",
                          return_value=self.empty_match) as patched_check:
            with patch.object(halogenases.HalogenaseResult, "get_best_match",
                              return_value=[self.test_trp_6_7_match]):
                record = DummyRecord(seq=test_fasta["ktzR"])
                positive_test = specific_analysis(record)

                patched_check.assert_called_once()
                assert positive_test and positive_test[0].position == 6

    # list of several matches
    @patch.object(subprocessing, "run_hmmsearch",
                  return_value=[FakeHit("start", "end", 500, "mibH")])
    @patch.object(secmet.Record, "get_cds_features_within_regions",
                  return_value=[DummyCDS(locus_tag="mibH", translation=test_fasta["mibH"])])
    def test_more_best_match(self, _patched_get_cds, run_hmmsearch, _patched_by_name):
        for value in run_hmmsearch.return_value:
            value.hsps = [FakeHSPHit("trp_5", "trp_5", bitscore=800)]
        with patch.object(halogenases, "check_for_halogenases",
                          return_value=self.empty_match) as patched_check:
            with patch.object(halogenases.HalogenaseResult, "get_best_match",
                              return_value=[self.test_trp_5_match, self.test_trp_6_7_match]):
                record = DummyRecord(seq=test_fasta["mibH"])
                positive_test = specific_analysis(record)
                assert len(positive_test) == 1
                assert patched_check.return_value.position is None
                assert positive_test and positive_test[0].position is None

    # no best match
    @patch.object(subprocessing, "run_hmmsearch", return_value=[FakeHit("start", "end", 500, "mibH")])
    @patch.object(secmet.Record, "get_cds_features_within_regions",
                  return_value=[DummyCDS(locus_tag="mibH", translation=test_fasta["mibH"])])
    def test_no_best_match(self, _patched_get_cds, run_hmmsearch, _patched_by_name):
        for value in run_hmmsearch.return_value:
            value.hsps = [FakeHSPHit("trp_5", "trp_5", bitscore=800)]
        with patch.object(halogenases, "check_for_halogenases", return_value=self.empty_match) as patched:
            with patch.object(halogenases.HalogenaseResult, "get_best_match",
                              return_value=[]):
                record = DummyRecord(seq=test_fasta["mibH"])
                positive_test = specific_analysis(record)
                assert patched.return_value.position is None
                assert positive_test and positive_test[0].position is None
