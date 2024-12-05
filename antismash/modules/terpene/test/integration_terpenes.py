# License: GNU Affero General Public License v3 or later
# A copy of GNU AGPL v3 should have been included in this software package in LICENSE.txt.

# for test files, silence irrelevant and noisy pylint warnings
# pylint: disable=use-implicit-booleaness-not-comparison,protected-access,missing-docstring

import tempfile
import unittest

import antismash
from antismash.common.secmet.test.helpers import DummyCDS, DummyRecord
from antismash.common.test.helpers import run_and_regenerate_results_for_module
from antismash.config import build_config, destroy_config, update_config
from antismash.modules import terpene


DATA = {
    "FQ762_RS00865": {
        "location": [173800, 174946],
        "translation": (
            "MPGLPPMVAPPAEQVDDPVTATDAANAVDAVLRGVLDERLRHCRAVDPLFARELADRLAALTARG"
            "GKRLRTAFAHCGWRAAGGSGDATAVLRTGAALELLQACALVHDDVMDGSVQRRGAPALHVDLARG"
            "HWAAGMHGSSESFGTSAAVLTGDLALAWADDLLTETALGTPHGPRLHGEWRAMRTEMVAGQYRDL"
            "HAQAARSSGVDEALAIATLKSALYTVARPLALGAVLAGAADGDALEALRAAGRCAGLAFQLRDDL"
            "LGAFGDPALTGKPADDDLRSRKLTYLLAVAVRLADAADDHLAAAALAPDADPKSEKAVRQVRSAL"
            "VRTGARDLVETKIGELTDMSLAHFDRCGARPAVRHEFAALIGRATGAVPRGTGEVV"
        ),
    },
    "FQ762_RS00875": {
        "location": [176500, 177496],
        "translation": (
            "MTRRELDAAGITDPLLRAAYTRCRRLNARHGKTYFLATRLLPLERRSAVHALYGFARWADDIVDD"
            "LDRHRTAGERDRALRRLESDLADGLRTGTGREPVVRAVADTADRYDIDHVLFADFLSSMRADLTV"
            "THYPTYADLQAYVHGSAAVIGLQMLPVLGTVTDRKEAAPHAAALGVAFQLTNFLRDVGEDLDRGR"
            "VYLPGSLLAAHGVDRPLLEWSRRTGRRDRRIRAALVAAEAMTRRVYRTAEPGIAMLDPRVRPCIR"
            "AAFTLYGGILDAIAEQDYTVLHRRAVVSRRRRAATAAEGVLCVAAARWRARTSPRVEAAGGPAAE"
            "WQGPVR"
        ),
    },
    "FQ762_RS00895": {
        "location": [180821, 182038],
        "translation": (
            "MTPRSAQDSDVLVIGGGAAGLSLAHRLTENGTAPAMTLVEPPDGPLRPAERTWCYWGAAADGLEE"
            "AVGASWSVLRLHGADGGSVTVDPAPFTYRMVRSADFERMVHGRLARTDGARLLRGTAESVRAVPA"
            "GTEVRCTLPGGRPLTLYARRVFDSRPLPELPPARTCLLQHFRGWFVHTRTDRFDPAVADLMDFRV"
            "PQPAHGLAFGYVLPLAPDRALVEYTEFSRAPLTTEAYESALGHYCRDILGLGELTVERTEQGVIP"
            "MTDARFPGRAGPAVYRIGTAGGATRPATGYTFAAVQRHSGAIAAALRDGHDRVPAPHGRRARAMD"
            "AVLLRALDTGRIDGPRFFTDLFRRVPAERLLRFLDGTTSLREEWGIGLRTPVRPMLRTAAEVPFL"
            "PRRSQPLARTGGNNR"
        ),
    },
}


FUNGAL_DATA = {
    "ATEG_03568": {
        "location": [2447, 5057],
        "translation": (
            "MDTISDVMKHCVPINYDDYDPLPADHFSTYPVFCSKATAEAIEASAEFTRKWKRACEVDGLHQDK"
            "LNFQACTTHLGHYNQWAYPDCLPERVSLNAVLSDCAFFWDGMSPVTGRSTNVMITSPVDVSDSIS"
            "AEKMNELTQDFGIAMLSELQSGRRIEPRFEINKMAVQVIRDFIAVDPFTGIGHLKEWKGHLDAQE"
            "KSTHNNMSWEKYVEHRVNESGGNWGISVGCWTNDIRISDEEKESVKYLTQLACAGGILGNDYYSF"
            "PKEFDEHHRSGTLDRLQNGVALLMREYGYTEEEAKEIIKKEVIIREKKWMDGFDAWSRQAGPETG"
            "EIRRYLVMTMALMSGSMFWMSHAGRYHRTDLATTAEDRATLIGKSRGALRVLAGYPPPKNLEGIV"
            "REPLASAVQDDNGHVQHKDAVADSSVRNGVHHAFKKRNSRNGKQNGTEGSKSTFTNGNYVQPAKL"
            "QQHGTSINSMAIYTAPFQEAAGDICDAPYSYIDSLPSKKNRNKLLDLLNDWLQVPPSSLKRIKNI"
            "VHMLHNSSLMLDDIEDASALRRGQPATHTFYGISQTINSANYVYVHAVHEVTRLYNPDADELRNL"
            "HRGQSLDLYWRHHARCPSMEEYIVMVDNKTGGLFRLMLRLMTAESSISRPLDTALCRLLTLTGRY"
            "YQIRDDYLNLASADYESKKGFCEDFDEGKFSLPLIHLLSHTRYPDRITSALFNRKPGTNLPYEMK"
            "RYILAEMEEVQTLAYSQDVLKYLHEELMHALDETENRLGANDGIRMMLLGMGPKLLLC"
        ),
    }
}


class TestTerpene(unittest.TestCase):
    def setUp(self):
        options = build_config(self.get_args(), isolated=True, modules=antismash.get_all_modules())
        self.options = update_config(options)
        terpene.prepare_data()

    def tearDown(self):
        destroy_config()

    def get_args(self):
        return ["--minimal", "--enable-terpene"]

    def test_full_pathway(self):
        features = []
        for name, data in DATA.items():
            features.append(DummyCDS(locus_tag=name, start=data["location"][0],
                                     end=data["location"][1], translation=data["translation"]))
        record = DummyRecord(seq="A" * features[-1].end, features=features)
        with tempfile.NamedTemporaryFile(suffix=".gbk") as temp:
            record.to_genbank(temp.name)
            results = run_and_regenerate_results_for_module(temp.name, terpene, self.options)
        assert results
        cluster_pred = results.cluster_predictions[1]
        assert len(cluster_pred.cds_predictions) == 3
        first_domains, second_domains, third_domains = cluster_pred.cds_predictions.values()

        assert len(first_domains) == 1
        domain = first_domains[0]
        assert domain.domain_type == "PT_FPPS_like"
        assert len(domain.reactions) == 1
        assert domain.reactions[0].substrates[0].name == "IPP"

        assert len(second_domains) == 1
        domain = second_domains[0]
        assert domain.domain_type == "PT_phytoene_like"
        assert len(domain.reactions) == 3
        assert domain.reactions[0].substrates[0].name == "FPP"
        assert domain.reactions[-1].products[0].name == "hydroxysqualene"

        assert len(third_domains) == 1
        domain = third_domains[0]
        assert domain.domain_type == "Lycopene_cycl"
        assert len(domain.reactions) == 1
        assert domain.reactions[0].substrates[0].name == "C40_carotenoid"

        assert cluster_pred.products[0].name == "C40_carotenoid"

    def test_fungal_pathway(self):
        features = []
        for name, data in FUNGAL_DATA.items():
            features.append(DummyCDS(locus_tag=name, start=data["location"][0],
                                     end=data["location"][1], translation=data["translation"]))
        record = DummyRecord(seq="A" * features[-1].end, features=features)
        with tempfile.NamedTemporaryFile(suffix=".gbk") as temp:
            record.to_genbank(temp.name)
            results = run_and_regenerate_results_for_module(temp.name, terpene, self.options)
        assert results
        cluster_pred = results.cluster_predictions[1]
        assert len(cluster_pred.cds_predictions) == 1
        first_domains = list(cluster_pred.cds_predictions.values())[0]

        assert len(first_domains) == 2
        domain1 = first_domains[0]
        assert domain1.domain_type == "T1TS"
        assert domain1.subtypes == ("T1TS_IV-V_a",)
        assert len(domain1.reactions) == 1
        assert domain1.reactions[0].substrates[0].name == "GFPP"

        domain2 = first_domains[1]
        assert domain2.domain_type == "PT_FPPS_like"
        assert domain2.subtypes == tuple()
        assert len(domain2.reactions) == 1
        assert domain2.reactions[0].substrates[0].name == "IPP"

        assert cluster_pred.products[0].name == "sesterterpene_C1-C15/C14-C18"
