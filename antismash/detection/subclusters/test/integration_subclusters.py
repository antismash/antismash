# License: GNU Affero General Public License v3 or later
# A copy of GNU AGPL v3 should have been included in this software package in LICENSE.txt.

# for test files, silence irrelevant and noisy pylint warnings
# pylint: disable=use-implicit-booleaness-not-comparison,protected-access,missing-docstring

import glob
import os
import tempfile
import unittest
from unittest.mock import patch

import antismash
from antismash.common import path
from antismash.common.secmet import Record
from antismash.common.secmet.locations import FeatureLocation
from antismash.common.secmet.test.helpers import DummyCDS, DummyRecord
from antismash.common.test.helpers import run_and_regenerate_results_for_module
from antismash.config import build_config, destroy_config, update_config
from antismash.detection import subclusters
from antismash.detection.subclusters import compounds, ruleset, signatures
from antismash.detection.subclusters.compounds import CompoundInfo


DATA = {
    # dhpg
    "AJAP_31990": {
        "location": [17866, 18657],
        "translation": (
            "MNGVRYEKKDHVAYVTMDRPEALNAMDRRMHAELAGIWDDVEADDDVRAVVLTGAGDRAFSVGQD"
            "LKERARLTDEGVEASTFGSSGQPGHPRLTDRFTLSKPVVARVHGYALGGGFELVLACDIVVASEE"
            "AVFGLPEVRLGLIAGAGGVFRLPRQLPQKVAMGYLLTGRRMDAATALRHGLVNEVVPFPELDRCV"
            "AEWTDDLVRAAPLSVRAIKEAALRSLDLPLEEAFKTSYPWEERRRTSADASEGARAFAEKRDPIW"
            "TGR"
        ),
    },
    # dhpg
    "AJAP_31995": {
        "location": [18654, 19868],
        "translation": (
            "MTALSPGLDLRALAGEAHRVDDEVRTLRAAFVEAHAEELYAELTDGRTRYLRIDELVRAAALAFP"
            "GLVPSEAQMAAERARPQAEKEGREIDQGIFLRGILRAPKAGPHLLDAMLRPTARARRLLPEFLET"
            "GVVQMEAVRLERRDGVAHLTLCRDDCLNAEDAQQVDDMETAVDLVLLDPSVRVGLLRGGVMSHPR"
            "YLGRRVFCAGINLKKLSSGDIPLVDFLLRRETGYIQKIFRGLLTDDSWHSRFTGKPWMAAVDSFA"
            "IGGGTQLLLVFDHVLAASDAYLSLPAAKEGIIPGVSNFRLSRIAGPRVARQVILGGRKLRADEPD"
            "ARSIVDEVVPPEEMDAAIDGALARLDGEAVAANRRMVNLSEEPPEEFRRYIAEFALQQALRIYGA"
            "DVIGKVDGFAVGSR"
        )
    },
    # dhpg
    "AJAP_32000": {
        "location": [18865, 20527],
        "translation": (
            "MVTRNEIKGELVLRFDGSRPLSAAAVEEIGAFCDRAEDQREPGPVTIHVTGAPPADWAKGLAVGL"
            "VSKWERVVRRFERLGRLTAAVASGECAGTALDLLLAADIRIAEPGTTLRLASAGGGTWPGMTVYR"
            "LTKQAGAAGIRRAVLLGAPIGTERALALDLIDEVSDDPAKTLAELDVVVDGAETAIRRQLIFEAG"
            "STTFEEALGSHLAAADRALRREAKS"
        )
    },
    # dhpg
    "AJAP_32005": {
        "location": [20521, 21621],
        "translation": (
            "MTAIIEPAEDLSVLTGLTEITRFAGVGTAVSESSYSQTELLEILDIEDPKIRSVFLNSAIDRRFL"
            "TLPPEGPGGARATEPQGDLLDKHKRIAVDMGCRALEACLKSAGATLSDLRHLCCVTSTGFLTPGL"
            "SALIIREMGIDPHCSRSDIVGMGCNAGLNALNVVAGWSAAHPGELGVVLCSEACSAAYALDGTMR"
            "TAVVNSLFGDGSAALALVSGDGRVPGPRVLKFASYIITDAVDAMRYDWDRDQDRFSFFLDPQIPY"
            "VVGAHAEIVIDRLLSGTGLRRSDIGHWLVHSGGKKVIDAVVVNLGLSRHDVRHTTGVLRDYGNLS"
            "SGSFLFSYERLADEDVARPGDYGVLMTMGPGSTIEMALIQW"
        )
    },
    # dhpg and hpg
    "AJAP_32060": {
        "location": [43933, 45294],
        "translation": (
            "MDNSRKRGVHRLFLTLRRTVEILVFMDSNGLSTHLNVETLHGSLTDPAISSMNLLNELIDEYPVA"
            "ISMAAGRPYEEFFDIRLIHTYLDAYCDHLRRDRKMDEAVVTRTLFQYGTTKGVIADLIAKNLAED"
            "ENIDAAPESVVVTVGAQEAMFLVLRTLRAGERDVLLAPAPTYVGLTGAALLTDTPVWPVRSNENG"
            "IDPDDLVLQLKRADEQGKRVRACYVTPNFANPTGTSMDLPSRHRLLDVAESNGILLLEDNAYGLF"
            "GSERLPSLKSLDRSGSVVYLGSFAKTGMPGARVGFAVADQRMADGGLFADQLSKLKGMLTVNTSP"
            "IAQAVIAGKLLLNDFSLTKANAREIAIYQRNLRLTLDALERGLGSCEGVSWNTPTGGFFVTVTVP"
            "FVVDDELLETAAREHGVLFTPMHHFYGGKGGFHQLRLSISLLTPELIEEGVARLAALIKPRLP"
        )
    },
    # hpg
    "AJAP_32035": {
        "location": [57846, 58916],
        "translation": (
            "MTHLCLDDLERAARTALPGEIWDFLAGGSGAEASLAANRTALDRIFVIPRMLRELTGGTTEAEVL"
            "GRRASLPVAVAPVAYQRLFHPEGELAAARAARDAGVPYTICTLSSVPLEEIAVVGGRPWFQLYWL"
            "RDEKRSLELVRRAEDAGCEAIVFTVDVPWMGRRLRDMRNGFALPESVTAANFDAGAAAHRRTEGL"
            "SAVADHTAREFAPATWESVEAVRAHTNLPVVLKGILAVEDAVRAVDAGATGIVVSNHGGRQLDGA"
            "VPGIAMLEEIADAVSGGCEVLLDGGIRTGGDVLKALALGASGVLIGRPFMWGLAADGQAGARQVL"
            "DLLAVELRNALGLAGCDSVSAARRLSTRYER"
        )
    },
    # hpg
    "AJAP_32040": {
        "location": [58913, 59965],
        "translation": (
            "MSSHENTQNFEIDYVEMYVANLEVAASGWLDKYDFSVTATDRSADHRSVTLRHSAIALVLTEPLS"
            "DRHPGATYLQTHGEGVADIALRTTDVAAAFEAAVKAGAKPLREPEKGVDSVLTATVSGFGDVVHT"
            "LIQSDVAEEAPRGKGGVELGVIDHFAVCLNAGDLGPTVAFYERALGFKQIFEEHIVVGAQAMNST"
            "VVQSTSGAVTLTLIEPDKTADPGQIDDFIKEHHGSGVQHIAFTSPDAVRAVKELSARGVEFLKTP"
            "DAYYDLLGERIELETHSLDDLRETKLLADEDHGGQLFQIFTASTHPRKTIFFEIIERQGAGTFGS"
            "SNIKALYEAVELERTGQSKLGPARR"
        )
    },
    # hpg
    "AJAP_32155": {
        "location": [79880, 80710],
        "translation": (
            "MTIEKALVVGTGLIGTSAALSLREKGVAVHLSDIDAQAVRVARELGAGREWTGEEVDLAVIAVPP"
            "QLVGERLAELQKRGAARAYTDVASVKVDPIADAERLGCDMTSYVPGHPLAGRERSGPAAARAGLF"
            "AGRPWALCPGPGTGAEALRLTRGLVALCGAEAVTVGAAEHDSAVALVSHAPHVAASAVAASLASG"
            "DDVALSLAGQGLRDVTRIAAGNPLLWRRILSGNAVPVAAVLDRIAADLAAAATALRAGDLDDLTE"
            "LLRRGVDGHGRIPAIG"
        )
    },
    # hpg_2 (far away)
    "AJAP_32060_2": {
        "location": [343933, 345294],
        "translation": (
            "MDNSRKRGVHRLFLTLRRTVEILVFMDSNGLSTHLNVETLHGSLTDPAISSMNLLNELIDEYPVA"
            "ISMAAGRPYEEFFDIRLIHTYLDAYCDHLRRDRKMDEAVVTRTLFQYGTTKGVIADLIAKNLAED"
            "ENIDAAPESVVVTVGAQEAMFLVLRTLRAGERDVLLAPAPTYVGLTGAALLTDTPVWPVRSNENG"
            "IDPDDLVLQLKRADEQGKRVRACYVTPNFANPTGTSMDLPSRHRLLDVAESNGILLLEDNAYGLF"
            "GSERLPSLKSLDRSGSVVYLGSFAKTGMPGARVGFAVADQRMADGGLFADQLSKLKGMLTVNTSP"
            "IAQAVIAGKLLLNDFSLTKANAREIAIYQRNLRLTLDALERGLGSCEGVSWNTPTGGFFVTVTVP"
            "FVVDDELLETAAREHGVLFTPMHHFYGGKGGFHQLRLSISLLTPELIEEGVARLAALIKPRLP"
        )
    },
    # hpg_2 (far away)
    "AJAP_32035_2": {
        "location": [357846, 358916],
        "translation": (
            "MTHLCLDDLERAARTALPGEIWDFLAGGSGAEASLAANRTALDRIFVIPRMLRELTGGTTEAEVL"
            "GRRASLPVAVAPVAYQRLFHPEGELAAARAARDAGVPYTICTLSSVPLEEIAVVGGRPWFQLYWL"
            "RDEKRSLELVRRAEDAGCEAIVFTVDVPWMGRRLRDMRNGFALPESVTAANFDAGAAAHRRTEGL"
            "SAVADHTAREFAPATWESVEAVRAHTNLPVVLKGILAVEDAVRAVDAGATGIVVSNHGGRQLDGA"
            "VPGIAMLEEIADAVSGGCEVLLDGGIRTGGDVLKALALGASGVLIGRPFMWGLAADGQAGARQVL"
            "DLLAVELRNALGLAGCDSVSAARRLSTRYER"
        )
    },
    # hpg_2 (far away)
    "AJAP_32040_2": {
        "location": [358913, 359965],
        "translation": (
            "MSSHENTQNFEIDYVEMYVANLEVAASGWLDKYDFSVTATDRSADHRSVTLRHSAIALVLTEPLS"
            "DRHPGATYLQTHGEGVADIALRTTDVAAAFEAAVKAGAKPLREPEKGVDSVLTATVSGFGDVVHT"
            "LIQSDVAEEAPRGKGGVELGVIDHFAVCLNAGDLGPTVAFYERALGFKQIFEEHIVVGAQAMNST"
            "VVQSTSGAVTLTLIEPDKTADPGQIDDFIKEHHGSGVQHIAFTSPDAVRAVKELSARGVEFLKTP"
            "DAYYDLLGERIELETHSLDDLRETKLLADEDHGGQLFQIFTASTHPRKTIFFEIIERQGAGTFGS"
            "SNIKALYEAVELERTGQSKLGPARR"
        )
    },
    # hpg_2 (far away)
    "AJAP_32155_2": {
        "location": [379880, 380710],
        "translation": (
            "MTIEKALVVGTGLIGTSAALSLREKGVAVHLSDIDAQAVRVARELGAGREWTGEEVDLAVIAVPP"
            "QLVGERLAELQKRGAARAYTDVASVKVDPIADAERLGCDMTSYVPGHPLAGRERSGPAAARAGLF"
            "AGRPWALCPGPGTGAEALRLTRGLVALCGAEAVTVGAAEHDSAVALVSHAPHVAASAVAASLASG"
            "DDVALSLAGQGLRDVTRIAAGNPLLWRRILSGNAVPVAAVLDRIAADLAAAATALRAGDLDDLTE"
            "LLRRGVDGHGRIPAIG"
        )
    }
}

# the core cluster detection finds a T3PKS protocluster around AJAP_32005
T3PKS_SPAN = FeatureLocation(521, 41621, 1)

# the subclusters detected in this record
DHPG_SPAN = FeatureLocation(17866, 45294, 1)  # overlaps with T3PKS_SPAN and HPG_SPAN
HPG_SPAN = FeatureLocation(43933, 80710, 1)  # overlaps with DHPG_SPAN
HPG_2_SPAN = FeatureLocation(343933, 380710, 1)  # overlaps with nothing

# fixed rules, signatures and compounds, so that changes to the real data files
# don't affect the tests
TEST_RULES = path.get_full_path(__file__, "data", "rules.txt")
TEST_DETAILS = path.get_full_path(__file__, "data", "hmmdetails.txt")
TEST_COMPOUNDS = {
    "TestDhpg": CompoundInfo(name="3,5-Dihydroxyphenylglycine (Dhpg)", smiles=None,
                             classification=["amino acid"]),
    "TestHpg": CompoundInfo(name="4-Hydroxyphenylglycine (Hpg)", smiles=None,
                            classification=["amino acid"]),
}


class TestSubclusters(unittest.TestCase):

    def setUp(self):
        # the aggregate built from the test signatures must not replace the real one
        temp_dir = tempfile.TemporaryDirectory()  # pylint: disable=consider-using-with
        self.addCleanup(temp_dir.cleanup)
        aggregate = os.path.join(temp_dir.name, "subcluster_seeds.hmm")
        test_signatures = {signature.name: signature
                           for signature in signatures._read_signatures(TEST_DETAILS)}
        patches = [
            patch.object(subclusters, "DETAILS_FILE", TEST_DETAILS),
            patch.object(subclusters, "AGGREGATE_HMM_FILE", aggregate),
            patch.object(ruleset, "AGGREGATE_HMM_FILE", aggregate),
            patch.dict(signatures._SIGNATURE_CACHE, test_signatures, clear=True),
            patch.object(ruleset, "_rule_files_for_strictness", return_value=[TEST_RULES]),
            patch.dict(ruleset._CACHE, clear=True),
            patch.dict(compounds._COMPOUND_CACHE, TEST_COMPOUNDS, clear=True),
        ]
        for patcher in patches:
            patcher.start()
            self.addCleanup(patcher.stop)
        subclusters.prepare_data()

    def tearDown(self):
        destroy_config()

    def get_args(self, mode="create"):
        return ["--minimal", "--subclusters", "--subclusters-region-mode", mode,
                "--enable-html"]

    def build_options(self, mode="create"):
        """ Replaces any active config with one for the given region mode """
        destroy_config()
        return update_config(build_config(self.get_args(mode), isolated=True,
                                          modules=antismash.get_all_modules()))

    def build_record(self):
        """ Builds the test record on its own, without any area from another module """
        features = []
        for name, data in DATA.items():
            features.append(DummyCDS(locus_tag=name, start=data["location"][0],
                                     end=data["location"][1], translation=data["translation"]))
        return DummyRecord(seq="A" * max(feature.end for feature in features), features=features)

    def run_with_mode(self, mode):
        """ Runs the whole pipeline in the given mode

            Returns this module's results along with the regions the pipeline
            formed, which have to be read back from the written output as the
            helper itself only returns the results.
        """
        options = self.build_options(mode)
        record = self.build_record()
        regions = []

        def capture_regions(output_dir):
            written = [name for name in glob.glob(os.path.join(output_dir, "*.gbk"))
                       if ".region" not in name]
            regions.extend(Record.from_genbank(written[0])[0].get_regions())

        with tempfile.NamedTemporaryFile(suffix=".gbk") as temp:
            record.to_genbank(temp.name)
            results = run_and_regenerate_results_for_module(temp.name, subclusters, options,
                                                            callback=capture_regions)
        return results, regions

    def check_detection(self, results):
        """ Checks the detection itself, which none of the region modes affect """
        assert len(results.rule_results.protoclusters) == 3
        protoclusters = results.rule_results.protoclusters
        assert [(proto.product, proto.location) for proto in protoclusters] == [
            ("TestDhpg", DHPG_SPAN),
            ("TestHpg", HPG_SPAN),
            ("TestHpg", HPG_2_SPAN),
        ]

        assert len(results.predictions) == 3
        dhpg, hpg, hpg_2 = results.predictions

        assert dhpg.rule.name == "TestDhpg"
        assert dhpg.compound.name == "3,5-Dihydroxyphenylglycine (Dhpg)"
        assert dhpg.location == DHPG_SPAN
        assert len(dhpg.domain_hits) == 5
        assert sorted((hit.cds_name, hit.domain_name) for hit in dhpg.domain_hits) == [
            ("AJAP_31990", "ECH_1"),
            ("AJAP_31995", "ECH_1"),
            ("AJAP_32000", "ECH_1"),
            ("AJAP_32005", "Chal_sti_synt_C"),
            ("AJAP_32060", "Aminotran_1_2"),
        ]
        aminotransferase = dhpg.domain_hits_by_cds["AJAP_32060"][0]
        assert aminotransferase.domain_accession.startswith("PF00155.")
        assert aminotransferase.domain_description == "Aminotransferase class I and II"

        assert hpg.rule.name == "TestHpg"
        assert hpg.compound.name == "4-Hydroxyphenylglycine (Hpg)"
        assert hpg.location == HPG_SPAN
        assert len(hpg.domain_hits) == 5
        assert sorted((hit.cds_name, hit.domain_name) for hit in hpg.domain_hits) == [
            ("AJAP_32035", "FMN_dh"),
            ("AJAP_32040", "Glyoxalase"),
            ("AJAP_32060", "Aminotran_1_2"),
            ("AJAP_32155", "PDH_C"),
            ("AJAP_32155", "PDH_N"),
        ]

        # same HPG subcluster, but further upstream
        assert hpg_2.location == HPG_2_SPAN

    def test_full_pathway_clip_mode(self):
        results, regions = self.run_with_mode("clip")
        assert results

        # detection itself, unaffected by the mode
        self.check_detection(results)

        # subregion formation: dhpg protocluster span clipped to T3PKS boundaries
        predicted_subregions = results.get_predicted_subregions()
        assert len(predicted_subregions) == 1
        assert predicted_subregions[0].label == "subclusters"
        assert predicted_subregions[0].tool == "subclusters"
        assert predicted_subregions[0].location == FeatureLocation(DHPG_SPAN.start,
                                                                   T3PKS_SPAN.end, 1)

        # region formation: only span of T3PKS protocluster
        assert len(regions) == 1
        assert regions[0].location == T3PKS_SPAN

        # predictions shown for the region: only dhpg overlaps it
        dhpg, _, _ = results.predictions
        assert results.get_predictions_for_region(regions[0]) == [dhpg]

    def test_full_pathway_extend_mode(self):
        results, regions = self.run_with_mode("extend")
        assert results

        # detection itself, unaffected by the mode
        self.check_detection(results)

        # subregion formation: full dhpg protocluster span
        predicted_subregions = results.get_predicted_subregions()
        assert len(predicted_subregions) == 1
        assert predicted_subregions[0].label == "subclusters"
        assert predicted_subregions[0].tool == "subclusters"
        assert predicted_subregions[0].location == DHPG_SPAN

        # region formation: spanning protoclusters for T3PKS and dhpg
        assert len(regions) == 1
        assert regions[0].location == FeatureLocation(T3PKS_SPAN.start, DHPG_SPAN.end, 1)

        # predictions shown for the region: dhpg, and hpg as it starts within dhpg
        dhpg, hpg, _ = results.predictions
        assert results.get_predictions_for_region(regions[0]) == [dhpg, hpg]

    def test_full_pathway_create_mode(self):
        results, regions = self.run_with_mode("create")
        assert results

        # detection itself, unaffected by the mode
        self.check_detection(results)

        # subregion formation: spanning protoclusters for dhpg and hpg
        predicted_subregions = results.get_predicted_subregions()
        assert len(predicted_subregions) == 2

        for subregion in predicted_subregions:
            assert subregion.label == "subclusters"
            assert subregion.tool == "subclusters"
        assert predicted_subregions[0].location == FeatureLocation(DHPG_SPAN.start, HPG_SPAN.end, 1)
        assert predicted_subregions[1].location == HPG_2_SPAN

        # region formation: spanning protoclusters for T3PKS, dhpg and hpg
        assert len(regions) == 2
        assert regions[0].location == FeatureLocation(T3PKS_SPAN.start, HPG_SPAN.end, 1)
        assert regions[1].location == HPG_2_SPAN

        # predictions shown for each region: the far away hpg only in the second
        dhpg, hpg, hpg_2 = results.predictions
        assert results.get_predictions_for_region(regions[0]) == [dhpg, hpg]
        assert results.get_predictions_for_region(regions[1]) == [hpg_2]
