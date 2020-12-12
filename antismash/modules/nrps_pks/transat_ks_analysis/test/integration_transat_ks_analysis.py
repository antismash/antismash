# License: GNU Affero General Public License v3 or later
# A copy of GNU AGPL v3 should have been included in this software package in LICENSE.txt.

# for test files, silence irrelevant and noisy pylint warnings
# pylint: disable=no-self-use,protected-access,missing-docstring

import unittest
from io import StringIO
from Bio import Phylo

from antismash.config import build_config, destroy_config
from antismash.common import path
from antismash.modules.nrps_pks.transat_ks_analysis.transat_ks_analysis import run_transpact_ks_analysis, \
    get_leaf2clade, _PPLACER_MASS_CUTOFF, _LEAF2CLADE_TBL, transpact_tree_prediction,  \
    get_transpact_clade


class TestTreePredictionMethods(unittest.TestCase):
    """ test all lose methods involved in predicting the transpact tree """

    def setUp(self):
        self.fun_clades, self.clade2ann = get_leaf2clade(_LEAF2CLADE_TBL)
        test_tree_path = path.get_full_path(__file__, "data", "test_pplacer_tree.tre")
        with open(test_tree_path, "r") as ffh:
            self.pplacer_tree = ffh.read()

        self.non_conserved_tree = "(((A_#0_M=1.0, etnangien_EtnE1_KS5_eDB)" \
                                  "(A_#0_M=1.0, oocydin_Ddad_OocN2_KS6_aDOHbLOH)))"

    def test_get_transpact_clade(self):
        """ test for transpact clade calling """
        newick_tree = Phylo.read(StringIO(self.pplacer_tree), 'newick')

        self.assertRaises(ValueError, get_transpact_clade, "not_present_query_name", newick_tree, self.fun_clades)

        non_conserved_newick_tree = Phylo.read(StringIO(self.non_conserved_tree), 'newick')
        clade_assignment = get_transpact_clade("A_#0_M=1.0", non_conserved_newick_tree, self.fun_clades)
        assert clade_assignment == "clade_not_conserved"

        clade_assignment = get_transpact_clade("nrpspksdomains_ctg1_11_PKS_KS.1_#0_M=1", newick_tree, self.fun_clades)
        assert clade_assignment == "Clade_30"

    def test_transpact_tree_prediction(self):
        """ test for transpact tree calling """
        prediction = transpact_tree_prediction(self.pplacer_tree, _PPLACER_MASS_CUTOFF, self.fun_clades, self.clade2ann)
        assert len(prediction.predictions) == 1
        assert prediction.predictions[0][0] == "glycine"
        assert prediction.predictions[0][1].clade == "Clade_30"
        assert prediction.predictions[0][1].mass_score == 1.0

        non_prediction = transpact_tree_prediction(self.non_conserved_tree, _PPLACER_MASS_CUTOFF, self.fun_clades,
                                                   self.clade2ann)
        assert len(non_prediction.predictions) == 1
        assert non_prediction.predictions[0][0] == "NA"
        assert non_prediction.predictions[0][1].clade == "clade_not_conserved"
        assert non_prediction.predictions[0][1].mass_score == 0.0

        self.assertRaises(ValueError, transpact_tree_prediction, "()", _PPLACER_MASS_CUTOFF, self.fun_clades,
                          self.clade2ann)


class TestRunTranspactKSAnalysis(unittest.TestCase):
    """ test cases """
    def setUp(self):
        build_config([])

        self.query_ks = {'nrpspksdomains_ctg1_11_PKS_KS.2': 'IAIIGLSGRYPEAENVNEFWDNLKTSKNCIREIPQERWDWKQYFDEEKGKAGKIYTKW'
                                                            'GGFLKEVDKFDPLFFQISPRDAEQMDPQERLFLEEAYKSIEDSGYTPETLCKSRKIGV'
                                                            'FVGVMNGTYAPQTRFFSMANRISYHLNFQGPSMAVDTACSASLTAIHLALESLYSGMS'
                                                            'ECAIAGGVNLILRPIHYMGLSAMTMLSSTNENRSFGDHADGFVDGEGVGAVVLKSLKA'
                                                            'AIADDDHIYGILKGSVVNAGGKTNGYTVPNPNAQYELVLEALKRSGIDARAVSCIEAH'
                                                            'GTGTVLGDPIEVSGLARAFEHYSKDRQFCTIGSVKSNIGHCESAAGIAGVTKVLLQMQ'
                                                            'QGQIVPSLHSKRLNPNIDFSNTPFVVQQELENWNRPVLKVDGQKKEYPRIAGISSFGA'
                                                            'GGANAHVVIEEY',
                         'nrpspksdomains_ctg1_11_PKS_KS.3': 'IAIVGLSGRYPQAPDLAAYWQNLRNGKDCITEVPPDRWDWKKYYTEDRNQGGRHYSKW'
                                                            'GGFIEDVDKFDPLFFNISPREAELIDPQERLFLEHVWMALEDAGYRGEDLQGEAGEYL'
                                                            'GGQVGVYAGVMYGEYQLFAAEESVRGNPLSVGGSYASIANRISYVLNLHGPSMSIDTM'
                                                            'CSSSLTALHLACQDLKLGRTNLGVAGGVNVTIHPNKYLMLSNGQFISSGGHCESFGKG'
                                                            'GDGYIPGEGVGVALLKRLSDAQRDGDHIYGVIKGSGLNHGGKTHGYSVPNPKAQQTVI'
                                                            'SRALKESGIRPGAISYIEAHGTGTKLGDPIEITGLSKAFGSEVQRGSCFIGSAKSNIG'
                                                            'HCESAAGIAGVTKVLLQMKYGQIVPSLHSQVLNPNIDFSQTPFVVQQELKEWFRPLLK'
                                                            'LNGTRKEYPRIAGVSSFGAGGANAHLVIEEY',
                         'nrpspksdomains_ctg1_12_PKS_KS.1': 'IAIVGMAGKYPGAKDLEEYWKNLVGGKNSIQEVPLSRWDVNQYYDPAPSKPGKVNCKW'
                                                            'LGLLDDIDCFDPLFFMISPAEAEEMDPQHRIFLEEGHKAFESAGYSSSGLSNKKCGVY'
                                                            'LGIMSNEYSLLAIEAQRSMSGTSNNYAIGAARMAYYLNLKGPAIPIDTACSSSLVGTH'
                                                            'LACQALHAGEIDMALVGGVTLYLTAESYIDMCASGMLSAEGQCKTFDDSADGFVPGEG'
                                                            'VGAIVLKRLSEAQRNDDCILGTIIGSGINQDGKTNGITAPSVKSQIELEREIYRKYKI'
                                                            'DPSSISYVETHGTGTKLGDPIELEALSTVFKEQTDRKNYCALGSVKSNIGHTSAAAGV'
                                                            'ASMHKALLSLKHQKLVPTLNFKKPNSHFNFDESPFYVNTELKRWESEGEKPRRVGVSS'
                                                            'FGFSGTNAHVVIEEY',
                         'nrpspksdomains_ctg1_12_PKS_KS.2': 'IAIIGLSGRYPQAYNVEEYWNHLKTGKDCITEIPRDRWSLEGFFHEDFEEAVAQGKSY'
                                                            'SKWGGFLEGFSEFDPLFFNISSHEAMSMDPQERLFLQTCWEVCEDGGYTREKIAIQHH'
                                                            'GKVGVFAGITKTGFDLYGPELWKKGEPVFPHTSFSSVANRVSYLLNLKGPSMPIDTMC'
                                                            'SSSLTAIHEACKHVDQGECELAIAGGVNLYLHPSTYRGLCAQRMLALDRQCRSFGKGG'
                                                            'NGFVPGEGVGAVLLKRLSQALFDGDHIYGVIRGTSANHGGKTNGYTVSNPKVQAELIR'
                                                            'ETLDKTGIDARTVSYIEAHGTATELGDPIEVAGLSQAFQQDTVDTGFCALGSAKSNVG'
                                                            'HLEAAAGMAGVTKILLQMKHKVLVPSLHAKELNPNINFSKTPFVVQQELGEWRSPFLE'
                                                            'RDGAGKEYPRIAGISSFGAGGANAHVIIEEY'}

    def tearDown(self):
        destroy_config()

    def test_run_transpact_ks_analysis(self):
        """ test the pipeline """
        results = run_transpact_ks_analysis(self.query_ks)

        assert 'nrpspksdomains_ctg1_11_PKS_KS.2' in results
        assert len(results['nrpspksdomains_ctg1_11_PKS_KS.2'].predictions) == 1
        assert results['nrpspksdomains_ctg1_11_PKS_KS.2'].predictions[0][0] == 'aMe_eDB'
        assert results['nrpspksdomains_ctg1_11_PKS_KS.2'].predictions[0][1].clade == 'Clade_47'
        assert results['nrpspksdomains_ctg1_11_PKS_KS.2'].predictions[0][1].mass_score == 1.0

        assert 'nrpspksdomains_ctg1_11_PKS_KS.3' in results
        assert len(results['nrpspksdomains_ctg1_11_PKS_KS.3'].predictions) == 1
        assert results['nrpspksdomains_ctg1_11_PKS_KS.3'].predictions[0][0] == 'eDB'
        assert results['nrpspksdomains_ctg1_11_PKS_KS.3'].predictions[0][1].clade == 'Clade_62'
        assert results['nrpspksdomains_ctg1_11_PKS_KS.3'].predictions[0][1].mass_score == 1.0

        assert 'nrpspksdomains_ctg1_12_PKS_KS.1' in results
        assert results['nrpspksdomains_ctg1_12_PKS_KS.1'].predictions
        assert results['nrpspksdomains_ctg1_12_PKS_KS.1'].predictions[0][0] == 'b_MeeDB'
        assert results['nrpspksdomains_ctg1_12_PKS_KS.1'].predictions[0][1].clade == 'Clade_82'
        assert results['nrpspksdomains_ctg1_12_PKS_KS.1'].predictions[0][1].mass_score == 1.0

        assert 'nrpspksdomains_ctg1_12_PKS_KS.2' in results
        assert len(results['nrpspksdomains_ctg1_12_PKS_KS.2'].predictions) == 1
        assert results['nrpspksdomains_ctg1_12_PKS_KS.2'].predictions[0][0] == 'clade_not_conserved'
        assert results['nrpspksdomains_ctg1_12_PKS_KS.2'].predictions[0][1].clade == 'NA'
        assert results['nrpspksdomains_ctg1_12_PKS_KS.2'].predictions[0][1].mass_score == 0.0
