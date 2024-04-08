# License: GNU Affero General Public License v3 or later
# A copy of GNU AGPL v3 should have been included in this software package in LICENSE.txt.

# for test files, silence irrelevant and noisy pylint warnings
# pylint: disable=use-implicit-booleaness-not-comparison,protected-access,missing-docstring

import json
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
from antismash.detection.genefunctions.halogenases.flavin_dependent.subgroups import (
    indolic,
    phenolic,
    pyrrolic
)

from antismash.detection.genefunctions.halogenases import halogenases_analysis
from antismash.detection.genefunctions.halogenases.halogenases import (
    HalogenaseHmmResult,
    FlavinDependentHalogenases, # add_potential_matches, get_best_match, finalize_enzyme
    Match # to_json, from_json
)

from antismash.detection.genefunctions.halogenases.flavin_dependent import substrate_analysis
from antismash.detection.genefunctions.halogenases.flavin_dependent.substrate_analysis import (
    run_halogenase_phmms,
    search_residues,
    categorize_on_substrate_level,
    check_for_halogenases,
    fdh_specific_analysis,
    retrieve_fdh_signature_residues,
    _gather_fdh_substrate_modules,
    _get_analysis_modules,
    _ANALYSIS_MODULES,
    _get_substrate_specific_profiles,
    check_for_match,
    search_conserved_motif #test
)

from antismash.detection.genefunctions.halogenases.flavin_dependent.substrate_analysis import (
    FDH_SUBGROUPS
)

""" KtzQ: Trp-7
    KtzR: Trp-6
    mibH: Trp-5
    bmp2: pyrrole
    ChlB4: orsellinic
    End30: hpg
    BhaA: tyrosine
"""

test_protein_translations = {
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
                APARHEFKLVGEQAAELARTLPTQYEYFAQLHRAR""",
        "bmp2": """MDQFKSYDVVIIGSGPAGSLCGIECRKKGLSVLCIEKDEFPRFHIGESLTGNAGQIIRDL
                GLADEMNAAGFPDKPGVNVIGSLSKNEFFIPILAPTWQVRRSDFDNMLKRRALEHGVEYQ
                QGLVKDVIKHEEKVVGAIYKADGVEHQVRSKVLVDASGQNTFLSRKGIAGKREIEFFSQQ
                IASFAHYKNVERDLPPFSTNTTILYSKQYHWSWIIPISPDTDSLGIVIPKDLYYKECKNP
                DDAIEWGMEHISPEIRRRFKNAERVGESQSMADFSYRIEPFVGDGWLCIGDAHRFLDPIF
                SYGVSFAMKEGIKAADAIKRAIDGNDWKTPFYEYRDWSNGGQQIAADLIRYFWIYPIFFG
                YQMQNPDLRDEVIRLLGGCCFDCEGWKAPTIFRNAIEEYDRKQMAG""",
        "ChlB4": """MQPDFDAAIVGGGPAGSAMASYLAEAGLSVAVFESEMFPRPHIGESLVPATMPVLDEIGV
                MPDIEAAGFPKKYGAAWTSAESRDVPHNGFTGLDHDFKAAEVMFVERDQPGVHRDYTFHV
                DRGKFDLILLKHAESRGAQVFQKTRVLKADFDTDPDLVTLNCRLGPRTLDFTTRMVIDAS
                GRQTMLGNQLKVKVPDPVFNQYAIHAWFEGLDRTAMALDPAKRDYIYVHFLPLEDTWMWQ
                IPITDTITSVGVVTQKHRFKAASADREKFFWDIVSSRKDIYDALQKAERIRPFKAEGDYS
                YAMRQICGDRFLLIGDAARFVDPIFSSGVSVALNSARLAAKDVIAAHRAGDFRKESFATY
                EEKLRRAVRNWYEFISVYYRLNILFTAFVQDPRYRIDVLKMLQGDFYDGEEPKALKAMRD
                LVTKVENDPEHLWHPYLGTLRAPSAAPTF""",
        "End30": """STLSTLVAMQGHSVLLLEKETFPRYQIGESLLPSTIHGICHLLGVTDELAAAGFPHKRGG
                TFRWGASPKPWNFSFSVSSKVSGPTSFAYQVERSKFDKILLDNAARKGVVVRQDRTVTDV
                VDDADGRARGLRYTDPDGTEHEVSARYVVDASGNTSRIHKRVGGSRTYSDFFKSLALFGY
                FENGKRMPAPYAGNILCVAFGSGWFWYIPLSSTLTSVGAVVRREDAAKVQGDPESALRGL
                IDECPMIKEYLADATRVTTGQYGQLRVRKDYSYHHTTFWRPGMVLVGDAACFVDPVFSSG
                VHLATYSALLAARSLNSVLAGRIDERRAFDEFEARYRREYGVFYEFLTSFYDMHVDEDSY
                FWTAKK""",
        "BhaA": """STVATLVAMQGHRVLLLEKEVFPRYQIGESLLPATVHGVCRMLGISDELANAGFPIKRGG
                TFRWGARPEPWTFHFGISAKMAGSTSHAYQVERARFDEMLLNNAKRKGVVVREGCAVTDV
                VEDGERVTGARYTDPDGTEREVSARFVIDASGNKSRLYTKVGGSRNYSEFFRSLALFGYF
                EGGKRLPEPVSGNILSVAFDSGWFWYIPLSDTLTSVGAVVRREDAEKIQGDREKALNTLI
                AECPLISEYLADATRVTTGRYGELRVRKDYSYQQETYWRPGMILVGDAACFVDPVFSSGV
                HLATYSALLAARSINSVLAGDLDEKTALNEFELRYRREYGVFYEFLVSFYQMNVNEESYF
                WQAKK""",
        "VatD": """MIEVCIIGFGFSAVPLIRELQRTGTEFKIISEESNSVWDALSQSNRLDFDLVSSYLTSFY
                SFDLVKDFVEDYYPTSKQFYEMHQRWRKVYENEIIRDRVTRIDNFKEHSVIFTKSGKTLN
                AKHVICSTGFSRAIHTHINDIDYSVSNKTFVFDTMGDSANLIISKLIPNNNKIIIRTNGF
                NARDKVVPGAGAIYTLDQLEGHNFRYMSHEHYGSVIYGLPIGSKKPILMGDQFPVTVRDD
                NYITSKSRPASGTIAIKYWPIDQYADKFGNNLEESISQGYLLNDIAMWLHTGKAIVVPKD
                TAINFEKKTITYAGIERAFHQYIKGDPEQPRLPKIMIDGNTPYEYQYRDNFMGVIPRTLN
                NVYFIGYTRPYTGGLANIIEMQGLFVHKMITQSEFHQKIHHNLDERIVAYNNHYYGTTKP
                RSADHLVYFGFYTDDLARLIGIDYKPSEINSIKDMVFYYAFPNNALKYRLKGEYAVDGVE
                DLIKKINEKYYDFIDVFAYLWGTSKMDSVELTEELEQFIRQYFNDMRHKEPYTKFLENYI
                QVYRRVKNTRVDETDDYEWSLMVKKASETRDRVLQEFKESGDYQLEENFRNFISNEIELI
                QSLMNYKILSVKDGQLKIVPEKIGGHSILENILCKITKILGFQGIIESVLGKNEMLPKIE
                AQDLQSLLSLTKPKEYELLYLKP""",
        "CtoA": """MEANPTAGTEVVVIGAGIVGVHSAIQFAKRGLKVVLIDNIVGQKKSFKVGESLLVFSNMF
                LRTISELDEFNQKCFPKHGVWFTYGMEGTTSFEEKAEWALESTLPQAMRDAFANKALLRA
                MADDVQIVRPEAEELMQQTARAHPNITFLDTAKVTNVVIAEGGGPHEVTWECKATQRTGV
                VRTTWLIDCSGRNRLLAKKLKHAAEDVELNDGFKTTAVWGQFSGIKDEMFGENWVNRTSD
                GARSKRDLNTLHLWGDGYWIWVIRLSEGRISVGATYDQRRPPAGAGYREKFWDIIRRYPL
                FDGMLSDDNMLEFHVFKNCQHITDTFVSEKRYGMIGDAASVIDAYYSQGVSLALVTSWHI
                TNIMERDLRERRLDKEYIARVNRHTRQDWHIMRNMVIEKYTSAMADGRFFVMTHLLDMII
                FVGAAFPRYLLVRWLVETQGSTARETPVMREMRRYLEENLYYSKIGSLAPEKVQKVQRGL
                QASLSERARWRVENGVKVTRLKAIVHAPSGLLKFWKLPLSGQREFEDISPKPVKQIPKWL
                AMTGEETNPRMLKMARPLMASTFFLMYGYDGLSTAVTKVRQRLERLPGAAATAETTAAGR
                RGEAPEPAMNGAAPVRNVLREVSA"""
        }

class TestPhenolic(unittest.TestCase):
    def setUp(self):
        self.test_tyrosine_match = Match("tyrosine-like_hpg_FDH", "flavin", "FHD",
                                         None, 1, None, "Tyr")
        self.test_hpg_match = Match("tyrosine-like_hpg_FDH", "flavin", "FDH",
                                    None, 1, None, "Hpg")
        self.test_cycline_orsellinic = Match("cycline_orsellinic_FDH", "flavin", "FDH",
                                             None, 1, None, "orsellinic")
        
        self.tyrosine_hmm_result = HalogenaseHmmResult(
            hit_id='tyrosine-like_hpg_FDH',
            bitscore=1000,
            query_id='tyrosine-like_hpg_FDH',
            enzyme_type="Flavin-dependent",
            profile=phenolic.SPECIFIC_PROFILES[1].path
            )
        self.hpg_hmm_result = HalogenaseHmmResult(
            hit_id='tyrosine-like_hpg_FDH',
            bitscore=600,
            query_id='tyrosine-like_hpg_FDH',
            enzyme_type="Flavin-dependent",
            profile=phenolic.SPECIFIC_PROFILES[1].path
            )
        self.cycline_orsellinic_hmm_result = HalogenaseHmmResult(
            hit_id='cycline_orsellinic_FDH',
            bitscore=600,
            query_id='cycline_orsellinic_FDH',
            enzyme_type="Flavin-dependent",
            profile=phenolic.SPECIFIC_PROFILES[0].path
            )

        self.tyr_empty_enzyme = FlavinDependentHalogenases("BhaA")
        self.tyr_enzyme_with_matches = FlavinDependentHalogenases("BhaA", potential_matches=self.test_tyrosine_match)
        
        self.hpg_empty_enzyme = FlavinDependentHalogenases("End30")
        self.hpg_enzyme_with_matches = FlavinDependentHalogenases("End30", potential_matches=self.test_hpg_match)
        
        self.orsellinic_empty_enzyme = FlavinDependentHalogenases("ChlB4")
        self.orsellinic_enzyme_with_matches = FlavinDependentHalogenases("ChlB4", potential_matches=self.test_cycline_orsellinic)
        
        phenolic_matches = [self.test_tyrosine_match, self.test_hpg_match, self.test_cycline_orsellinic]
        
    def test_check_for_fdh(self):
        with patch.object(substrate_analysis, "get_residues",
                          return_value={"cycline_orsellinic_FDH": phenolic.OTHER_PHENOLIC_SIGNATURE_RESIDUES}):
            categorize_on_substrate_level(DummyCDS(), self.orsellinic_empty_enzyme, self.cycline_orsellinic_hmm_result)
            assert self.orsellinic_empty_enzyme.potential_matches[0].profile == "cycline_orsellinic_FDH"
            assert self.orsellinic_empty_enzyme.potential_matches[0].confidence == 1

        with patch.object(substrate_analysis, "get_residues",
                          return_value={"tyrosine-like_hpg_FDH": phenolic.HPG_SIGNATURE_RESIDUES}):
            categorize_on_substrate_level(DummyCDS(), self.hpg_empty_enzyme, self.hpg_hmm_result)
            assert self.hpg_empty_enzyme.potential_matches[0].profile == "tyrosine-like_hpg_FDH"
            
        with patch.object(substrate_analysis, "get_residues",
                          return_value={"tyrosine-like_hpg_FDH": phenolic.TYROSINE_LIKE_SIGNATURE_RESIDUES}):
            categorize_on_substrate_level(DummyCDS(), self.tyr_empty_enzyme, self.tyrosine_hmm_result)
            assert self.tyr_empty_enzyme.potential_matches[0].profile == "tyrosine-like_hpg_FDH"

        with patch.object(substrate_analysis, "get_residues",
                          return_value={"tyrosine-like_hpg_FDH": None}):
            assert categorize_on_substrate_level(DummyCDS(), self.tyr_empty_enzyme,
                                 self.tyrosine_hmm_result) is None

        with patch.object(substrate_analysis, "get_residues",
                          return_value={"tyrosine-like_hpg_FDH": phenolic.TYROSINE_LIKE_SIGNATURE_RESIDUES}):
            substrate = FlavinDependentHalogenases("tyrosine-like_hpg_FDH", phenolic.TYROSINE_LIKE_SIGNATURE_RESIDUES)
            with patch.object(substrate_analysis, "TailoringEnzymes",
                              return_value=substrate):
                low_quality_hit = HalogenaseHmmResult("wrong_name", 400, "tyrosine-like_hpg_FDH", "foo", "tyrosine-like_hpg_FDH")
                assert not categorize_on_substrate_level(DummyCDS(), self.tyr_empty_enzyme, low_quality_hit)

    def test_get_signatures(self):
        assert phenolic.get_signatures() == [phenolic.TYROSINE_LIKE_SIGNATURE,
                                             phenolic.HPG_SIGNATURE,
                                             phenolic.OTHER_PHENOLIC_SIGNATURE]
        
    def test_get_residues(self):
        cds = DummyCDS(locus_tag="BhaA", translation=test_protein_translations["BhaA"])
        residues = substrate_analysis.retrieve_fdh_signature_residues(cds.translation, self.hpg_hmm_result,
                                                   [phenolic.HPG_SIGNATURE, phenolic.TYROSINE_LIKE_SIGNATURE],
                                                   enzyme_substrates=["Tyr", "Hpg"])
        assert isinstance(residues, dict)
        
class TestPyrrolic(unittest.TestCase):
    def setUp(self):
        self.pyrrole_empty_enzyme = FlavinDependentHalogenases("bmp2")
        self.test_pyrrole_match = Match("pyrrole_FDH", "flavin", "FDH",
                                  "mono/di", 1, None, "pyrrole")
        
        self.pyrrole_hmm_result = HalogenaseHmmResult(
            hit_id='pyrrole_FDH',
            bitscore=600,
            query_id='pyrrole_FDH',
            enzyme_type="Flavin-dependent",
            profile=pyrrolic.SPECIFIC_PROFILES[0].path
        )

        self.pyrrole_enzyme_with_matches = FlavinDependentHalogenases("bmp2",
                                                            potential_matches=[self.test_pyrrole_match])
        
    def test_check_for_fdh(self):
        with patch.object(substrate_analysis, "get_residues",
                          return_value=pyrrolic.PYRROLE_SIGNATURE_RESIDUES):
            categorize_on_substrate_level(DummyCDS(), self.pyrrole_empty_enzyme, self.pyrrole_hmm_result)
            assert self.pyrrole_empty_enzyme.potential_matches[0].profile == "pyrrole_FDH"
            assert self.pyrrole_empty_enzyme.potential_matches[0].confidence == 1

        with patch.object(substrate_analysis, "get_residues",
                          return_value=None):
            assert categorize_on_substrate_level(DummyCDS(), self.pyrrole_empty_enzyme,
                                 self.pyrrole_hmm_result) is None

        with patch.object(substrate_analysis, "get_residues",
                          return_value=pyrrolic.PYRROLE_SIGNATURE_RESIDUES):
            substrate = FlavinDependentHalogenases("pyrrole_FDH", pyrrolic.PYRROLE_SIGNATURE_RESIDUES)
            with patch.object(substrate_analysis, "TailoringEnzymes",
                              return_value=substrate):
                low_quality_hit = HalogenaseHmmResult("wrong_name", 400,
                                                      "pyrrole_FDH", "foo", "pyrrole_FDH")
                assert not categorize_on_substrate_level(DummyCDS(), self.pyrrole_empty_enzyme, low_quality_hit)
                
    def test_get_residues(self):
        cds = DummyCDS(locus_tag="bmp2", translation=test_protein_translations["bmp2"])
        residues = substrate_analysis.retrieve_fdh_signature_residues(cds.translation, self.pyrrole_hmm_result,
                                                   [pyrrolic.PYRROLE_SIGNATURE],
                                                   enzyme_substrates=["mono_di", "tetra", "unconv_mono_di"])
        assert isinstance(residues, dict)


class TestIndolic(unittest.TestCase):
    def setUp(self):
        self.test_trp_5_match = Match("trp_5_FDH", "flavin", "FDH",
                                      5, 0.5, None, "tryptophan")
        self.test_trp_6_7_match = Match("trp_6_7_FDH", "flavin", "FDH",
                                        6, 1, "", "tryptophan")
        
        self.trp_5_hmm_result = HalogenaseHmmResult(
            hit_id='trp_5_FDH',
            bitscore=1000,
            query_id='trp_5_FDH',
            enzyme_type='Flavin-dependent',
            profile=indolic.SPECIFIC_PROFILES[0].path,
        )
        self.trp_6_7_hmm_result = HalogenaseHmmResult(
            hit_id='trp_6_7_FDH',
            bitscore=1000,
            query_id='trp_6_7_FDH',
            enzyme_type='Flavin-dependent',
            profile=indolic.SPECIFIC_PROFILES[1].path,
        )

        tryptophan_matches = [self.test_trp_5_match, self.test_trp_6_7_match]
        
        self.trp_enzyme_with_matches = FlavinDependentHalogenases("ktzR", potential_matches=tryptophan_matches)
        self.trp_empty_enzyme = FlavinDependentHalogenases("ktzQ")

    def test_get_best_match(self):
        assert not self.trp_empty_enzyme.get_best_match()

        positive_test_best_match = self.trp_enzyme_with_matches.get_best_match()
        assert len(positive_test_best_match) == 1 and isinstance(positive_test_best_match[0], Match)

        self.trp_enzyme_with_matches.add_potential_matches(self.test_trp_6_7_match)
        assert len(self.trp_enzyme_with_matches.potential_matches) == 3

        multiple_matches = self.trp_enzyme_with_matches.get_best_match()
        assert len(multiple_matches) == 2 and isinstance(positive_test_best_match[0], Match)

        one_potential_match = FlavinDependentHalogenases("test_enzyme",
                                               potential_matches=[self.test_trp_5_match])
        multiple_matches = one_potential_match.get_best_match()
        assert len(multiple_matches) == 1 and isinstance(positive_test_best_match[0], Match)
        
# halogenases analysis
    def test_conversion_methods(self):
        converted_to_json = self.trp_enzyme_with_matches.to_json()
        converted_to_json = json.loads(json.dumps(converted_to_json))
        assert isinstance(converted_to_json, dict)

        converted_from_json = FlavinDependentHalogenases.from_json(converted_to_json)
        assert isinstance(converted_from_json, FlavinDependentHalogenases)
        assert converted_to_json == converted_from_json.to_json()

    @patch.object(subprocessing, "run_hmmsearch",
                  return_value=[FakeHit("start", "end", 1000, "foo")])
    def test_run_halogenase_phmms(self, run_hmmsearch):
        for value in run_hmmsearch.return_value:
            value.hsps = [FakeHSPHit("foo", "foo", bitscore=250)]

        negative_test_halogenase_hmms_by_id = run_halogenase_phmms("", [])
        assert not negative_test_halogenase_hmms_by_id

        for value in run_hmmsearch.return_value:
            value.hsps = [FakeHSPHit("foo", "foo", bitscore=1000)]

        positive_test_halogenase_hmms_by_id = run_halogenase_phmms("", [])
        for hit in positive_test_halogenase_hmms_by_id["foo"]:
            assert isinstance(hit, HalogenaseHmmResult)

    def test_search_signature_residues(self):
        positions = indolic.TRP_6_SIGNATURE

        with patch.object(subprocessing.hmmpfam, "run_hmmpfam2") as run_hmmpfam2:
            run_hmmpfam2.return_value = []
            signature_residues = search_residues(test_protein_translations["ktzR"],
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

            signature_residues = search_residues(list(test_protein_translations.values())[0],
                                                           positions, self.trp_6_7_hmm_result)
            assert signature_residues is None

            # checking if hit_id == query_id it runs til the reference position extraction
            for result in run_hmmpfam2.return_value:
                result.hsps = [FakeHSPHit("trp_6_7_FDH", hit_id="trp_6_7_FDH")]
                for hit in result.hsps:
                    hit_query = DummyFeature()
                    hit_profile = DummyFeature()
                    hit_query.seq = test_protein_translations["ktzR"]
                    hit_profile.seq = "xfdhcgkbjlnkml"
                    hit.aln = [hit_profile, hit_query]

            with patch.object(utils, "extract_by_reference_positions", return_value="dummy"):
                signature_residues = search_residues(test_protein_translations["ktzR"],
                                                               positions,
                                                               self.trp_6_7_hmm_result)
                assert signature_residues == "dummy"

    def test_false_check_for_match(self):
        false_test = FlavinDependentHalogenases("fake_name")
        false_test_residues = "FAKESEQUENCE"
        FDH_SUBGROUPS["trp_5_FDH"].update_match("trp_5_FDH", false_test_residues,
                                                false_test, self.trp_6_7_hmm_result)
        assert false_test.family == "" and false_test.cofactor == ""

    def test_check_for_fdh(self):
        with patch.object(substrate_analysis, "get_residues",
                          return_value={"trp_5_FDH": indolic.TRP_5_SIGNATURE_RESIDUES}):
            categorize_on_substrate_level(DummyCDS(), self.trp_empty_enzyme, self.trp_5_hmm_result)
            assert self.trp_empty_enzyme.potential_matches[0].profile == "trp_5_FDH"
            assert self.trp_empty_enzyme.potential_matches[0].confidence == 1
            assert self.trp_empty_enzyme.potential_matches[0].position == 5

        with patch.object(substrate_analysis, "get_residues",
                          return_value={"trp_5_FDH": indolic.TRP_5_SIGNATURE_RESIDUES}):
            low_quality_hit = HalogenaseHmmResult("trp_5_FDH", 380, "trp_5_FDH", "foo", "trp_5_FDH")
            categorize_on_substrate_level(DummyCDS(), self.trp_empty_enzyme, low_quality_hit)
            assert self.trp_empty_enzyme.potential_matches[1].profile == "trp_5_FDH"
            assert self.trp_empty_enzyme.potential_matches[1].confidence == 0.5

        with patch.object(substrate_analysis, "get_residues",
                          return_value={"trp_6_7_FDH": indolic.TRP_6_SIGNATURE_RESIDUES}):
            categorize_on_substrate_level(DummyCDS(), self.trp_empty_enzyme, self.trp_6_7_hmm_result)
            assert self.trp_empty_enzyme.potential_matches[2].profile == "trp_6_7_FDH"
            assert self.trp_empty_enzyme.potential_matches[2].position == 6
            
        with patch.object(substrate_analysis, "get_residues", return_value={"trp_6_7_FDH": "VISIL"}):
            categorize_on_substrate_level(DummyCDS(), self.trp_empty_enzyme, self.trp_6_7_hmm_result)
            assert self.trp_empty_enzyme.potential_matches[3].profile == "trp_6_7_FDH"
            assert self.trp_empty_enzyme.potential_matches[3].position == 7

        with patch.object(substrate_analysis, "get_residues",
                          return_value={"trp_6_7_FDH": None}):
            assert categorize_on_substrate_level(DummyCDS(), self.trp_empty_enzyme,
                                 self.trp_6_7_hmm_result) is None

        with patch.object(substrate_analysis, "get_residues",
                          return_value={"trp_5_FDH": indolic.TRP_5_SIGNATURE_RESIDUES}):
            substrate = FlavinDependentHalogenases("trp_5_FDH", indolic.TRP_5_SIGNATURE_RESIDUES)
            with patch.object(substrate_analysis, "TailoringEnzymes",
                              return_value=substrate):
                low_quality_hit = HalogenaseHmmResult("wrong_name", 400, "trp_5_FDH", "foo", "trp_5_FDH")
                assert not categorize_on_substrate_level(DummyCDS(), self.trp_empty_enzyme, low_quality_hit)

    def test_check_for_halogenases(self):
        negative_checked_halogenases = check_for_halogenases(DummyCDS(), [])
        assert negative_checked_halogenases is None

        mock_cds_feature = DummyCDS(locus_tag="foo", translation="")
        positive_checked_halogenase = check_for_halogenases(mock_cds_feature, [self.trp_5_hmm_result])
        assert isinstance(positive_checked_halogenase, FlavinDependentHalogenases)
        assert positive_checked_halogenase.cds_name == "foo"

    def test_get_signatures(self):
        assert indolic.get_signatures() == [indolic.TRP_5_SIGNATURE,
                                            indolic.TRP_6_SIGNATURE]

@patch.object(secmet.Record, "get_cds_by_name", return_value="ktzR")
class TestSpecificAnalysis(unittest.TestCase):
    def setUp(self):
        self.record = DummyRecord()
        self.cds = DummyCDS(locus_tag="ktzR", translation=test_protein_translations["ktzR"])
        self.test_trp_5_match = Match("trp_5_FDH", "flavin", "FDH", 5, 0, "", "tryptophan")
        self.test_trp_6_7_match = Match("trp_6_7_FDH", "flavin", "FDH", 6, 1, "", "tryptophan")
        self.empty_match = FlavinDependentHalogenases("test_enzyme")
        
    # one best match
    @patch.object(subprocessing, "run_hmmsearch",
                  return_value=[FakeHit("start", "end", 900, "ktzR")])
    @patch.object(secmet.Record, "get_cds_features_within_regions",
                  return_value=[DummyCDS(locus_tag="ktzR", translation=test_protein_translations["ktzR"])])
    def test_one_best_match(self, _patched_get_cds, run_hmmsearch, _patched_by_name):
        for value in run_hmmsearch.return_value:
            value.hsps = [FakeHSPHit("trp_6_7_FDH", "trp_6_7_FDH", bitscore=1500)]
        with patch.object(substrate_analysis, "check_for_halogenases",
                          return_value=self.empty_match) as patched_check:
            with patch.object(FlavinDependentHalogenases, "get_best_match",
                              return_value=[self.test_trp_6_7_match]):
                record = DummyRecord(seq=test_protein_translations["ktzR"])
                positive_test = fdh_specific_analysis(record)
                patched_check.assert_called_once()
                assert positive_test[0].substrates == "tryptophan"
                assert positive_test and positive_test[0].target_positions == 6

    # list of several matches
    @patch.object(subprocessing, "run_hmmsearch",
                  return_value=[FakeHit("start", "end", 500, "mibH")])
    @patch.object(secmet.Record, "get_cds_features_within_regions",
                  return_value=[DummyCDS(locus_tag="mibH", translation=test_protein_translations["mibH"])])
    def test_more_best_match(self, _patched_get_cds, run_hmmsearch, _patched_by_name):
        for value in run_hmmsearch.return_value:
            value.hsps = [FakeHSPHit("trp_5_FDH", "trp_5_FDH", bitscore=800)]
        with patch.object(substrate_analysis, "check_for_halogenases",
                          return_value=self.empty_match) as patched_check:
            with patch.object(FlavinDependentHalogenases, "get_best_match",
                              return_value=[self.test_trp_5_match, self.test_trp_6_7_match]):
                record = DummyRecord(seq=test_protein_translations["mibH"])
                positive_test = fdh_specific_analysis(record)
                assert len(positive_test) == 1
                assert patched_check.return_value.target_positions is None
                assert positive_test and positive_test[0].target_positions is None

    # no best match
    @patch.object(subprocessing, "run_hmmsearch", return_value=[FakeHit("start", "end", 500, "mibH")])
    @patch.object(secmet.Record, "get_cds_features_within_regions",
                  return_value=[DummyCDS(locus_tag="mibH", translation=test_protein_translations["mibH"])])
    def test_no_best_match(self, _patched_get_cds, run_hmmsearch, _patched_by_name):
        for value in run_hmmsearch.return_value:
            value.hsps = [FakeHSPHit("trp_5_FDH", "trp_5_FDH", bitscore=800)]
        with patch.object(substrate_analysis, "check_for_halogenases", return_value=self.empty_match) as patched:
            with patch.object(FlavinDependentHalogenases, "get_best_match",
                              return_value=[]):
                record = DummyRecord(seq=test_protein_translations["mibH"])
                positive_test = fdh_specific_analysis(record)
                assert patched.return_value.target_positions is None
                assert positive_test and positive_test[0].target_positions is None


class TestGeneralEnzymes(unittest.TestCase):
    def setUp(self):
        self.record = DummyRecord()
        self.general_cds = DummyCDS(locus_tag="CtoA", translation=test_protein_translations["CtoA"])
        self.general_match = Match("all_general_FDH", "flavin", "FDH", None, None, None)
        self.general_empty_match = FlavinDependentHalogenases("CtoA")
        
        self.unconventional_cds = DummyCDS(locus_tag="VatD", translation=test_protein_translations["VatD"])
        self.unconventional_match = Match("unconventional_FDH", "flavin", "FDH", None, None, None)
        self.unconventional_empty_match = FlavinDependentHalogenases("VatD")

   
    @patch.object(secmet.Record, "get_cds_by_name",
                  return_value=DummyCDS(locus_tag="VatD", translation=test_protein_translations["VatD"]))
    @patch.object(secmet.Record, "get_cds_features_within_regions",
                  return_value=[DummyCDS(locus_tag="VatD", translation=test_protein_translations["VatD"])])
    def test_unconventional_fdh_specific_analysis(self, _patched_get_cds, _patched_get_cds_by_name):
        record = DummyRecord(seq=test_protein_translations["VatD"])
        positive_test = fdh_specific_analysis(record)
        assert positive_test[0].substrates is None
        assert not positive_test[0].consensus_residues
        

    
    @patch.object(secmet.Record, "get_cds_by_name",
                  return_value=DummyCDS(locus_tag="CtoA", translation=test_protein_translations["CtoA"]))
    @patch.object(secmet.Record, "get_cds_features_within_regions",
                  return_value=[DummyCDS(locus_tag="CtoA", translation=test_protein_translations["CtoA"])])
    def test_conventional_fdh_specific_analysis(self, _patched_get_cds, _patched_get_cds_by_name):
        record = DummyRecord(seq=test_protein_translations["CtoA"])
        positive_test = fdh_specific_analysis(record)
        assert positive_test[0].substrates is None
        assert positive_test[0].consensus_residues == {'W.W.I.': 'WIWVIR'}
        
        with patch.object(substrate_analysis, "search_signature_residues",
                      return_value=None):
            record = DummyRecord(seq=test_protein_translations["CtoA"])
            positive_test = fdh_specific_analysis(record)
            assert not positive_test[0].consensus_residues

    @patch.object(secmet.Record, "get_cds_by_name",
                  return_value=DummyCDS(locus_tag="CtoA", translation=test_protein_translations["CtoA"]))
    @patch.object(secmet.Record, "get_cds_features_within_regions",
                  return_value=[DummyCDS(locus_tag="CtoA", translation=test_protein_translations["CtoA"])])
    def test_specific_analysis(self, _patched_get_cds, _patched_get_cds_by_name):
        record = DummyRecord(seq=test_protein_translations["CtoA"])
        categorized_halogenase = halogenases_analysis.specific_analysis(record)
        assert categorized_halogenase is not None