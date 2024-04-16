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
from antismash.common.subprocessing import hmmscan

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
    FlavinDependentHalogenases,
    Match
)

from antismash.detection.genefunctions.halogenases.flavin_dependent import substrate_analysis
from antismash.detection.genefunctions.halogenases.flavin_dependent.substrate_analysis import (
    run_halogenase_phmms,
    search_residues,
    categorize_on_substrate_level,
    fdh_specific_analysis
)

from antismash.detection.genefunctions.halogenases.flavin_dependent.substrate_analysis import (
    FDH_SUBGROUPS
)

#
""" KtzQ: Trp-7
    KtzR: Trp-6
    mibH: Trp-5
    bmp2: Pyrrole
    ChlB4: Orsellinic-like
    End30: Hpg
    BhaA: Tyrosine
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
        "CtoA": """MEANPTAGTEVVVIGAGIVGVHSAIQFAKRGLKVVLIDNIVGQKKSFKVGESLLVFSNMFLRTISELDEFNQKCFPKHGVWFTYGMEGTTSFEEKAEWALESTLPQAMRDAFANKALLRAMADDVQIVRPEAEELMQQTARAHPNITFLDTAKVTNVVIAEGGGPHEVTWECKATQRTGVVRTTWLIDCSGRNRLLAKKLKHAAEDVELNDGFKTTAVWGQFSGIKDEMFGENWVNRTSDGARSKRDLNTLHLWGDGYWIWVIRLSEGRISVGATYDQRRPPAGAGYREKFWDIIRRYPLFDGMLSDDNMLEFHVFKNCQHITDTFVSEKRYGMIGDAASVIDAYYSQGVSLALVTSWHITNIMERDLRERRLDKEYIARVNRHTRQDWHIMRNMVIEKYTSAMADGRFFVMTHLLDMIIFVGAAFPRYLLVRWLVETQGSTARETPVMREMRRYLEENLYYSKIGSLAPEKVQKVQRGLQASLSERARWRVENGVKVTRLKAIVHAPSGLLKFWKLPLSGQREFEDISPKPVKQIPKWLAMTGEETNPRMLKMARPLMASTFFLMYGYDGLSTAVTKVRQRLERLPGAAATAETTAAGRRGEAPEPAMNGAAPVRNVLREVSA
        """
        }

test_protein_translations = {key: "".join(val.split()) for key, val in test_protein_translations.items()}

class PhenolicBase(unittest.TestCase):
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
        self.less_confident_tyrosine_hmm_result = HalogenaseHmmResult(
            hit_id='tyrosine-like_hpg_FDH',
            bitscore=310,
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

        # tyrosine
        self.tyr_empty_enzyme = FlavinDependentHalogenases("BhaA", "flavin", "FDH")
        self.tyr_enzyme_with_matches = FlavinDependentHalogenases("BhaA", "flavin", "FDH",
                                                                  potential_matches=self.test_tyrosine_match)

        # Hpg
        self.hpg_empty_enzyme = FlavinDependentHalogenases("End30", "flavin", "FDH")
        self.hpg_enzyme_with_matches = FlavinDependentHalogenases("End30", "flavin", "FDH",
                                                                  potential_matches=self.test_hpg_match)

        # other phenolic-substrate halogenase (orsellinic-like)
        self.orsellinic_empty_enzyme = FlavinDependentHalogenases("ChlB4", "flavin", "FDH")
        self.orsellinic_enzyme_with_matches = FlavinDependentHalogenases("ChlB4", "flavin", "FDH",
                                                                         potential_matches=self.test_cycline_orsellinic)

        self.tyrosine_match_not_hpg = {"Tyr": "GFQRLGDAGLSGVPSYGADPSGLYW",
                                        "Hpg": "VALAMI"}

        phenolic_matches = [self.test_tyrosine_match, self.test_hpg_match,
                            self.test_cycline_orsellinic]

class PyrrolicBase(unittest.TestCase):
    def setUp(self):
        self.pyrrole_empty_enzyme = FlavinDependentHalogenases("bmp2", "flavin", "FDH")
        self.test_pyrrole_match = Match("pyrrole_FDH", "flavin", "FDH",
                                  "tetra", 1, None, "pyrrole")

        self.pyrrole_hmm_result = HalogenaseHmmResult(
            hit_id='pyrrole_FDH',
            bitscore=1000,
            query_id='pyrrole_FDH',
            enzyme_type="FDH",
            profile=pyrrolic.SPECIFIC_PROFILES[0].path
        )

        self.negative_pyrrole_hmm_result = HalogenaseHmmResult(
            hit_id='pyrrole_FDH',
            bitscore=1,
            query_id='pyrrole_FDH',
            enzyme_type="FDH",
            profile=pyrrolic.SPECIFIC_PROFILES[0].path
        )

        # tetrabrominating halogenase
        self.pyrrole_enzyme_with_matches = FlavinDependentHalogenases("bmp2", "flavin", "FDH",
                                                                      potential_matches=[self.test_pyrrole_match])

class IndolicBase(unittest.TestCase):
    def setUp(self):
        self.test_trp_5_match = Match("trp_5_FDH", "flavin", "FDH", 1, None,
                                      "tryptophan", 5, "mono")
        self.test_trp_6_7_match = Match("trp_6_7_FDH", "flavin", "FDH", 1, "",
                                      "tryptophan", 6, "mono")

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

        tryptophan_single_matches = [self.test_trp_6_7_match]
        tryptophan_matches = [self.test_trp_5_match, self.test_trp_6_7_match]

        # Trp-5 halogenase
        self.trp_5_enzyme_with_matches = FlavinDependentHalogenases("mibH", "flavin", "FDH",
                                                                    potential_matches=tryptophan_matches)
        # Trp-6 halogenase
        self.trp_enzyme_with_matches = FlavinDependentHalogenases("ktzR", "flavin", "FDH",
                                                                  potential_matches=tryptophan_single_matches)
        # Trp-7 halogenase
        self.trp_empty_enzyme = FlavinDependentHalogenases("ktzQ", "flavin", "FDH")

class TestPhenolic(PhenolicBase):
    def test_categorize_on_substrate_level(self):
        # Other phenolic (orsellinic-like) substrate halogenases
        with patch.object(substrate_analysis, "retrieve_fdh_signature_residues",
                    return_value=""):
            below_cutoff_hit = HalogenaseHmmResult("wrong_name", 400, "wrong_name", "foo", "wrong_name")
            assert not categorize_on_substrate_level(DummyCDS(), FlavinDependentHalogenases("ChlB4", "flavin", "FDH"),
                                                     [below_cutoff_hit])

        with patch.object(substrate_analysis, "retrieve_fdh_signature_residues",
                          return_value={"cycline_orsellinic_FDH": phenolic.OTHER_PHENOLIC_SIGNATURE_RESIDUES}):
            categorize_on_substrate_level(DummyCDS(), self.orsellinic_empty_enzyme, [self.cycline_orsellinic_hmm_result])
            assert self.orsellinic_empty_enzyme.potential_matches[0].profile == "cycline_orsellinic_FDH"
            assert self.orsellinic_empty_enzyme.potential_matches[0].confidence == 1
            assert self.orsellinic_empty_enzyme.potential_matches[0].substrates == "cycline_orsellinic-like"
            assert self.orsellinic_empty_enzyme.potential_matches[0].target_positions == [6, 8]

        with patch.object(substrate_analysis, "retrieve_fdh_signature_residues",
                        return_value={"cycline_orsellinic_FDH": ""}):
            assert not categorize_on_substrate_level(DummyCDS(), self.orsellinic_empty_enzyme,
                                                     [self.cycline_orsellinic_hmm_result])

        # Tyr and Hpg halogenases
        with patch.object(substrate_analysis, "retrieve_fdh_signature_residues",
                          return_value=phenolic.TYR_HPG_SIGNATURE_RESIDUES):
            categorize_on_substrate_level(DummyCDS(), self.hpg_empty_enzyme, [self.hpg_hmm_result])
            assert self.hpg_empty_enzyme.potential_matches[0].profile == "tyrosine-like_hpg_FDH"
            assert self.hpg_empty_enzyme.potential_matches[0].confidence == 1
            assert self.hpg_empty_enzyme.potential_matches[0].substrates == "Hpg"
            assert self.hpg_empty_enzyme.potential_matches[0].target_positions == [6,8]

        with patch.object(substrate_analysis, "retrieve_fdh_signature_residues",
                          return_value=self.tyrosine_match_not_hpg):
            assert categorize_on_substrate_level(DummyCDS(), self.tyr_empty_enzyme,
                                 [self.less_confident_tyrosine_hmm_result])
            assert self.tyr_empty_enzyme.potential_matches[0].profile == "tyrosine-like_hpg_FDH"
            assert self.tyr_empty_enzyme.potential_matches[0].confidence == 0.5
            assert self.tyr_empty_enzyme.potential_matches[0].substrates == "Tyr"
            assert self.tyr_empty_enzyme.potential_matches[0].target_positions == [6,8]

        with patch.object(substrate_analysis, "retrieve_fdh_signature_residues",
                          return_value={None}):
            assert categorize_on_substrate_level(DummyCDS(), self.tyr_empty_enzyme,
                                 self.tyrosine_hmm_result) is None

    def test_negative_search_for_match(self):
        assert not phenolic.search_for_match(phenolic.TYR_HPG_SIGNATURE_RESIDUES,
                                         self.tyr_empty_enzyme,
                                         self.less_confident_tyrosine_hmm_result,
                                         [6,8],
                                         1000, expected_residues=phenolic.TYR_HPG_SIGNATURE_RESIDUES)

        assert not phenolic.search_for_match(phenolic.TYR_HPG_SIGNATURE_RESIDUES,
                                         self.tyr_empty_enzyme,
                                         self.less_confident_tyrosine_hmm_result,
                                         [6,8],
                                         [1000, 500], expected_residues=phenolic.TYR_HPG_SIGNATURE_RESIDUES)

        assert not phenolic.search_for_match(phenolic.TYR_HPG_SIGNATURE_RESIDUES,
                                         self.tyr_empty_enzyme,
                                         self.less_confident_tyrosine_hmm_result,
                                         [6,8],
                                         "false_cutoff", expected_residues=phenolic.TYR_HPG_SIGNATURE_RESIDUES)

    def test_retrieve_fdh_signature_residues(self):
        cds = DummyCDS(locus_tag="BhaA", translation=test_protein_translations["BhaA"])
        assert " " not in cds.translation
        residues = substrate_analysis.retrieve_fdh_signature_residues(cds.translation, self.hpg_hmm_result,
                                                   [phenolic.HPG_SIGNATURE, phenolic.TYROSINE_LIKE_SIGNATURE],
                                                   enzyme_substrates=["Tyr", "Hpg"])
        # returns None for each substrate, but works, when running antismash
        self.fail(residues)
        assert isinstance(residues, dict)

class TestPyrrolic(PyrrolicBase):
    def test_get_consensus_signature(self):
        assert pyrrolic.get_consensus_signature(DummyCDS(translation=test_protein_translations["bmp2"]),
                                                self.pyrrole_hmm_result)

    def test_check_for_fdh(self):
        with patch.object(pyrrolic, "get_consensus_signature",
                          return_value={"pyrrole_FDH": pyrrolic.PYRROLE_SIGNATURE_RESIDUES}) as patched:
            assert not categorize_on_substrate_level(DummyCDS(), self.pyrrole_empty_enzyme,
                                                     [self.negative_pyrrole_hmm_result])
            patched.assert_called_once()

        with patch.object(pyrrolic, "get_consensus_signature",
                          return_value={"pyrrole_FDH": None}) as patched:
            assert not categorize_on_substrate_level(DummyCDS(), self.pyrrole_empty_enzyme,
                                                     [self.pyrrole_hmm_result])
            patched.assert_called_once()

        # conventional mono/dihalogenating pyrrole-halogenase
        with patch.object(pyrrolic, "get_consensus_signature",
                    return_value={"pyrrole_FDH": "DRSVFW"}) as patched:
            assert categorize_on_substrate_level(DummyCDS(), self.pyrrole_empty_enzyme,
                                                 [self.pyrrole_hmm_result])
            assert self.pyrrole_empty_enzyme.potential_matches[0].profile == "pyrrole_FDH"
            assert self.pyrrole_empty_enzyme.potential_matches[0].confidence == 1
            assert self.pyrrole_empty_enzyme.potential_matches[0].substrates == "pyrrole"
            assert self.pyrrole_empty_enzyme.potential_matches[0].number_of_decorations == "mono_di"
            patched.assert_called_once()

        # unconventional mono/dihalogenating pyrrole-halogenase
        with patch.object(pyrrolic, "get_consensus_signature",
                    return_value={"pyrrole_FDH": "YRRNFN"}) as patched:
            assert categorize_on_substrate_level(DummyCDS(), self.pyrrole_empty_enzyme,
                                                 [self.pyrrole_hmm_result])
            assert self.pyrrole_empty_enzyme.potential_matches[1].profile == "pyrrole_FDH"
            assert self.pyrrole_empty_enzyme.potential_matches[1].confidence == 1
            assert self.pyrrole_empty_enzyme.potential_matches[1].substrates == "pyrrole"
            assert self.pyrrole_empty_enzyme.potential_matches[1].number_of_decorations == "unconv_mono_di"
            patched.assert_called_once()

        # tetrahalogenating pyrrole-halogenase
        with patch.object(pyrrolic, "get_consensus_signature",
                    return_value={"pyrrole_FDH": "RRYFFA"}) as patched:
            assert categorize_on_substrate_level(DummyCDS(), self.pyrrole_empty_enzyme,
                                                 [self.pyrrole_hmm_result])
            assert self.pyrrole_empty_enzyme.potential_matches[2].profile == "pyrrole_FDH"
            assert self.pyrrole_empty_enzyme.potential_matches[2].confidence == 1
            assert self.pyrrole_empty_enzyme.potential_matches[2].substrates == "pyrrole"
            assert self.pyrrole_empty_enzyme.potential_matches[2].number_of_decorations == "tetra"
            patched.assert_called_once()

    def test_negative_search_for_match(self):
        assert not pyrrolic.search_for_match(pyrrolic.PYRROLE_SIGNATURE_RESIDUES,
                                         self.pyrrole_empty_enzyme,
                                         self.pyrrole_hmm_result,
                                         1000, expected_residues={"mono_di":"",
                                                             "unconv_mono_di":"",
                                                             "tetra":""})

class TestIndolic(IndolicBase):
    def test_get_best_match(self):
        assert not self.trp_empty_enzyme.get_best_match()

        positive_test_best_match = self.trp_enzyme_with_matches.get_best_match()

        assert len(positive_test_best_match) == 1 and isinstance(positive_test_best_match[0], Match)

        self.trp_enzyme_with_matches.add_potential_matches(self.test_trp_6_7_match)
        assert len(self.trp_enzyme_with_matches.potential_matches) == 2

        multiple_matches = self.trp_enzyme_with_matches.get_best_match()
        assert len(multiple_matches) == 2 and isinstance(positive_test_best_match[0], Match)

        one_potential_match = FlavinDependentHalogenases("test_enzyme", "flavin", "FDH",
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
        target_positions = indolic.TRP_6_SIGNATURE

        with patch.object(subprocessing.hmmpfam, "run_hmmpfam2") as run_hmmpfam2:
            run_hmmpfam2.return_value = []
            signature_residues = search_residues(test_protein_translations["ktzR"],
                                                           target_positions, self.trp_6_7_hmm_result)
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
                                                           target_positions, self.trp_6_7_hmm_result)
            assert signature_residues is None

            # checking if hit_id == query_id it runs til the reference target_positions extraction
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
                                                               target_positions,
                                                               self.trp_6_7_hmm_result)
                assert signature_residues == "dummy"

    def test_false_check_for_match(self):
        false_test = FlavinDependentHalogenases("fake_name", "flavin", "FDH")
        false_test_residues = "FAKESEQUENCE"
        FDH_SUBGROUPS["trp_5_FDH"].update_match("trp_5_FDH", false_test_residues,
                                                false_test, self.trp_6_7_hmm_result)
        assert false_test.substrates == None and false_test.target_positions == None and not false_test.potential_matches

    def test_check_for_fdh(self):
        with patch.object(indolic, "get_consensus_signature",
                          return_value={"trp_5_FDH": indolic.TRP_5_SIGNATURE_RESIDUES}):
            categorize_on_substrate_level(DummyCDS(), self.trp_empty_enzyme, [self.trp_5_hmm_result])
            assert self.trp_empty_enzyme.potential_matches[0].profile == "trp_5_FDH"
            assert self.trp_empty_enzyme.potential_matches[0].confidence == 1
            assert self.trp_empty_enzyme.potential_matches[0].substrates == "tryptophan"
            assert self.trp_empty_enzyme.potential_matches[0].number_of_decorations == "mono"
            assert self.trp_empty_enzyme.potential_matches[0].target_positions == 5

        with patch.object(indolic, "get_consensus_signature",
                          return_value={"trp_5_FDH": indolic.TRP_5_SIGNATURE_RESIDUES}):
            low_quality_hit = HalogenaseHmmResult("trp_5_FDH", 380, "trp_5_FDH", "foo", "trp_5_FDH")
            categorize_on_substrate_level(DummyCDS(), self.trp_empty_enzyme, [low_quality_hit])
            assert self.trp_empty_enzyme.potential_matches[1].profile == "trp_5_FDH"
            assert self.trp_empty_enzyme.potential_matches[1].confidence == 0.5
            assert self.trp_empty_enzyme.potential_matches[1].substrates == "tryptophan"
            assert self.trp_empty_enzyme.potential_matches[1].number_of_decorations == "mono"
            assert self.trp_empty_enzyme.potential_matches[1].target_positions == 5

        with patch.object(indolic, "get_consensus_signature",
                          return_value={"trp_6_7_FDH": indolic.TRP_6_SIGNATURE_RESIDUES}):
            categorize_on_substrate_level(DummyCDS(), self.trp_empty_enzyme, [self.trp_6_7_hmm_result])
            assert self.trp_empty_enzyme.potential_matches[2].profile == "trp_6_7_FDH"
            assert self.trp_empty_enzyme.potential_matches[2].confidence == 1
            assert self.trp_empty_enzyme.potential_matches[2].target_positions == 6
            assert self.trp_empty_enzyme.potential_matches[2].substrates == "tryptophan"
            assert self.trp_empty_enzyme.potential_matches[2].number_of_decorations == "mono"

        with patch.object(substrate_analysis, "retrieve_fdh_signature_residues",
                          return_value={"trp_6_7_FDH": "VISIL"}):
            categorize_on_substrate_level(DummyCDS(), self.trp_empty_enzyme, [self.trp_6_7_hmm_result])
            assert self.trp_empty_enzyme.potential_matches[3].profile == "trp_6_7_FDH"
            assert self.trp_empty_enzyme.potential_matches[3].confidence == 1
            assert self.trp_empty_enzyme.potential_matches[3].target_positions == 7
            assert self.trp_empty_enzyme.potential_matches[3].substrates == "tryptophan"
            assert self.trp_empty_enzyme.potential_matches[3].number_of_decorations == "mono"

        with patch.object(substrate_analysis, "retrieve_fdh_signature_residues",
                          return_value={"trp_6_7_FDH": None}):
            assert categorize_on_substrate_level(DummyCDS(), self.trp_empty_enzyme,
                                 self.trp_6_7_hmm_result) is None

        with patch.object(substrate_analysis, "retrieve_fdh_signature_residues",
                          return_value=""):
                low_quality_hit = HalogenaseHmmResult("wrong_name", 400, "wrong_name",
                                                      "foo", "wrong_name")
                assert not categorize_on_substrate_level(DummyCDS(),
                                                         FlavinDependentHalogenases("mibH", "flavin", "FDH"),
                                                         [low_quality_hit])

    def test_categorize_on_substrate_level(self):
        negative_checked_halogenases = categorize_on_substrate_level(DummyCDS(), self.trp_empty_enzyme, [])
        assert negative_checked_halogenases is None

        mock_cds_feature = DummyCDS(locus_tag="mibH", translation=test_protein_translations["mibH"])
        positive_checked_halogenase = categorize_on_substrate_level(mock_cds_feature, FlavinDependentHalogenases("mibH", "flavin", "FDH"),
                                                                    [self.trp_5_hmm_result])
        assert isinstance(positive_checked_halogenase, FlavinDependentHalogenases)
        assert positive_checked_halogenase.cds_name == "mibH"
        assert positive_checked_halogenase.potential_matches == [Match(profile='trp_5_FDH', cofactor='flavin', family='FDH',
                                                                       confidence=1.0, consensus_residues='VSILIREPGLPRGVPRAVLPGEA',
                                                                       substrates='tryptophan', target_positions=5, number_of_decorations='mono')]

class TestSpecificAnalysis(IndolicBase):
    def setUp(self):
        super().setUp()
        self.record = DummyRecord()

    def get_result(self, name):
        for value in run_hmmsearch.return_value:
            value.hsps = [FakeHSPHit("trp_6_7_FDH", name, bitscore=1500)]
        with patch.object(substrate_analysis, "categorize_on_substrate_level",
                          return_value=self.trp_enzyme_with_matches) as patched_check:
            with patch.object(FlavinDependentHalogenases, "get_best_match",
                              return_value=[self.test_trp_6_7_match]):
                cds = DummyCDS(locus_tag=name,
                               translation=test_protein_translations[name])
                record = self.record
                record.add_cds_feature(cds)
                with patch.object(secmet.Record, "get_cds_features_within_regions", return_value=[cds]):
                    return fdh_specific_analysis(record)

    def TEST_TEMP(self):
        for value in run_hmmsearch.return_value:
            value.hsps = [FakeHSPHit("trp_6_7_FDH", "ktzR", bitscore=1500)]
        with patch.object(substrate_analysis, "categorize_on_substrate_level",
                          return_value=self.trp_enzyme_with_matches) as patched_check:
            with patch.object(FlavinDependentHalogenases, "get_best_match",
                              return_value=[self.test_trp_6_7_match]):
                cds = DummyCDS(locus_tag="ktzR",
                                         translation=test_protein_translations["ktzR"])
                record = self.record
                record.add_cds_feature(cds)
                with patch.object(secmet.Record, "get_cds_features_within_regions", return_value=[cds]):
                    positive_test = self.get_result("ktzR")
            patched_check.assert_called_once()
        assert positive_test[0].family == "FDH"
        assert positive_test[0].cds_name == "ktzR"
        assert positive_test[0].cofactor == "flavin"
        assert positive_test[0].substrates == "tryptophan"
        assert positive_test[0].number_of_decorations == "mono"
        assert positive_test[0].potential_matches
        assert positive_test and positive_test[0].target_positions == 6

    # one best match
    @patch.object(subprocessing, "run_hmmsearch",
                  return_value=[FakeHit("start", "end", 900, "ktzR")])
    @patch.object(secmet.Record, "get_cds_features_within_regions",
                  return_value=[DummyCDS(locus_tag="ktzR",
                                         translation=test_protein_translations["ktzR"])])
    def test_one_best_match(self, _patched_get_cds, run_hmmsearch):
        for value in run_hmmsearch.return_value:
            value.hsps = [FakeHSPHit("trp_6_7_FDH", "ktzR", bitscore=1500)]
        with patch.object(substrate_analysis, "categorize_on_substrate_level",
                          return_value=self.trp_enzyme_with_matches) as patched_check:
            with patch.object(FlavinDependentHalogenases, "get_best_match",
                              return_value=[self.test_trp_6_7_match]):
                cds = DummyCDS(locus_tag="ktzR",
                                         translation=test_protein_translations["ktzR"])
                record = self.record
                record.add_cds_feature(cds)
                with patch.object(secmet.Record, "get_cds_features_within_regions", return_value=[cds]):
                    positive_test = fdh_specific_analysis(record)
            patched_check.assert_called_once()
        assert positive_test[0].family == "FDH"
        assert positive_test[0].cds_name == "ktzR"
        assert positive_test[0].cofactor == "flavin"
        assert positive_test[0].substrates == "tryptophan"
        assert positive_test[0].number_of_decorations == "mono"
        assert positive_test[0].potential_matches
        assert positive_test and positive_test[0].target_positions == 6

    # list of several matches, no best match
    @patch.object(subprocessing, "run_hmmsearch",
                  return_value=[FakeHit("start", "end", 500, "mibH")])
    @patch.object(secmet.Record, "get_cds_features_within_regions",
                  return_value=[DummyCDS(locus_tag="mibH", translation=test_protein_translations["mibH"])])
    def test_more_best_match(self, _patched_get_cds, run_hmmsearch):
        for value in run_hmmsearch.return_value:
            value.hsps = [FakeHSPHit("trp_5_FDH", "mibH", bitscore=800)]
        with patch.object(substrate_analysis, "categorize_on_substrate_level",
                          return_value=self.trp_5_enzyme_with_matches) as patched_check:
            with patch.object(FlavinDependentHalogenases, "get_best_match",
                              return_value=[self.test_trp_5_match, self.test_trp_6_7_match]) as _patched:
                record = self.record
                positive_test = fdh_specific_analysis(record)
                _patched.assert_called_once()
                assert len(positive_test) == 1
                assert positive_test[0].potential_matches
                assert positive_test[0].consensus_residues is None
                assert positive_test[0].number_of_decorations is None
                assert positive_test[0].target_positions is None

class TestGeneralEnzymes(unittest.TestCase):
    def setUp(self):
        self.record = DummyRecord()
        self.general_cds = DummyCDS(locus_tag="CtoA",
                                    translation=test_protein_translations["CtoA"])
        self.general_match = Match("all_general_FDH", "flavin", "FDH", None, None, None)
        self.general_empty_match = FlavinDependentHalogenases("CtoA", "flavin", "FDH")

        self.unconventional_cds = DummyCDS(locus_tag="VatD",
                                           translation=test_protein_translations["VatD"])
        self.unconventional_match = Match("unconventional_FDH", "flavin", "FDH", None, None, None)
        self.unconventional_empty_match = FlavinDependentHalogenases("VatD", "flavin", "FDH")

    @patch.object(secmet.Record, "get_cds_by_name",
                  return_value=DummyCDS(locus_tag="VatD",
                                        translation=test_protein_translations["VatD"]))
    @patch.object(secmet.Record, "get_cds_features_within_regions",
                  return_value=[DummyCDS(locus_tag="VatD",
                                         translation=test_protein_translations["VatD"])])
    def test_unconventional_fdh_specific_analysis(self, _patched_get_cds, _patched_get_cds_by_name):
        record = DummyRecord()
        positive_test = fdh_specific_analysis(record)
        self.fail(positive_test)
        assert isinstance(positive_test[0], FlavinDependentHalogenases)
        assert positive_test[0].substrates is None
        assert not positive_test[0].consensus_residues
        assert not positive_test[0].potential_matches

    @patch.object(secmet.Record, "get_cds_by_name",
                  return_value=DummyCDS(locus_tag="CtoA",
                                        translation=test_protein_translations["CtoA"]))
    @patch.object(secmet.Record, "get_cds_features_within_regions",
                  return_value=[DummyCDS(locus_tag="CtoA",
                                         translation=test_protein_translations["CtoA"])])
    def test_conventional_fdh_specific_analysis(self, _patched_get_cds, _patched_get_cds_by_name):
        record = DummyRecord(seq=test_protein_translations["CtoA"])
        positive_test = fdh_specific_analysis(record)
        assert positive_test[0].substrates is None
        assert not positive_test[0].potential_matches
        assert positive_test[0].consensus_residues == {'W.W.I.': 'WIWVIR'}

        with patch.object(substrate_analysis, "retrieve_fdh_signature_residues",
                      return_value=None):
            record = DummyRecord(seq=test_protein_translations["CtoA"])
            positive_test = fdh_specific_analysis(record)
            assert not positive_test[0].substrates and not positive_test[0].target_positions

    @patch.object(secmet.Record, "get_cds_by_name",
                  return_value=DummyCDS(locus_tag="CtoA",
                                        translation=test_protein_translations["CtoA"]))
    @patch.object(secmet.Record, "get_cds_features_within_regions",
                  return_value=[DummyCDS(locus_tag="CtoA",
                                         translation=test_protein_translations["CtoA"])])
    def test_specific_analysis(self, _patched_get_cds, _patched_get_cds_by_name):
        record = DummyRecord(seq=test_protein_translations["CtoA"])
        categorized_halogenase = halogenases_analysis.specific_analysis(record)


    @patch.object(secmet.Record, "get_cds_by_name",
                  return_value=DummyCDS(locus_tag="VatD",
                                        translation=test_protein_translations["VatD"]))
    @patch.object(secmet.Record, "get_cds_features_within_regions",
                  return_value=[DummyCDS(locus_tag="VatD",
                                         translation=test_protein_translations["VatD"])])
    @patch.object(substrate_analysis, "search_residues",
                  return_value=None)
    def test_negative_get_conserved_motifs(self, _patched_search_residues,
                                           _patched_get_cds, _patched_get_cds_by_name):
        record = DummyRecord(seq=test_protein_translations["VatD"])
        categorized_halogenase = halogenases_analysis.specific_analysis(record)
        assert categorized_halogenase[0].consensus_residues == {}


    @patch.object(secmet.Record, "get_cds_by_name",
                  return_value=DummyCDS(locus_tag="VatD",
                                        translation=test_protein_translations["VatD"]))
    @patch.object(secmet.Record, "get_cds_features_within_regions",
                  return_value=[DummyCDS(locus_tag="VatD",
                                         translation=test_protein_translations["VatD"])])
    @patch.object(hmmscan, "run_hmmscan",
                  return_value=[FakeHSPHit("any", "any")])
    def test_negative_enzyme_hits_fasta(self, _patched_run_hmmscan, _patched_get_cds, _patched_get_cds_by_name):
        for value in _patched_run_hmmscan.return_value:
            value.id = "test_id"
            value.hits = "test_for_pass"
        record = DummyRecord(seq=test_protein_translations["VatD"])
        assert not fdh_specific_analysis(record)
