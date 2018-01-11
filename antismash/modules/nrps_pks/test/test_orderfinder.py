# License: GNU Affero General Public License v3 or later
# A copy of GNU AGPL v3 should have been included in this software package in LICENSE.txt.

# for test files, silence irrelevant and noisy pylint warnings
# pylint: disable=no-self-use,protected-access,missing-docstring

import os
import unittest

from antismash.common import path
from antismash.common.test.helpers import DummyCDS
from antismash.detection.nrps_pks_domains import domain_identification
from antismash.modules.nrps_pks import orderfinder


class TestCTerminalExtract(unittest.TestCase):
    def setUp(self):
        # 110 long, since only the last 100 will be used
        self.seqs = {"STAUR_3982": ("FLEFTRQRGFISEEFGREHDSELMKTYLPTLRKDLVLLESYSYAEEAPLDMPLTV"
                                    "FASTRDRIIPSTQLESWGELTREKPSIHLFEGDHFFARDAGGPLLALIREKLGLG"),
                     "STAUR_3984": ("LARVLRMEASRIDRLRALGELGLDSLMSLELRNRLEASLGLKLSVTLLFTYPNLA"
                                    "GLAEYLHGELLPAAAREQPAAQSQTHAAPSQIAEQVEQLSKDELLAFFDKSFGIA"),
                     "STAUR_3983": ("TNMGLDSLMSLELRNRLEATLGLKLSATLLFTYPNLAALADHLLGKLSSVDEAPA"
                                    "KTAPTAAAPPPPPTLKPQAALPAELDQLGKDELLSLFDESLTESLKRTRMTRTSR"),
                     "STAUR_3985": ("PSKIDRLRALGELGLDSLMSLELRNRLEAALGMKLSATLLFTYPNLASLAQHVVG"
                                    "RMEFPSEATVAPITASPGAVEGQAERLAEVEQMSDDEAEQLLLASLESLSTELLK"),
                     "STAUR_3972": ("SEAALRGSAAGVAYTASKHALIGFTKNTAFMYGAKGVRVNIVAPGPVRTSISGAS"
                                    "RSDHGWSRIAPVMNVLAVPVAESATLAGHILWLMSDEAENINGAVLPSDGGWSTF")
                     }

        self.data_dir = path.get_full_path(os.path.dirname(__file__), "data", "terminals")

        self.features = [DummyCDS(1, 200, translation=seq, locus_tag=seq_id) for seq_id, seq in self.seqs.items()]
        self.features_by_id = {feature.locus_tag: feature for feature in self.features}

    def test_c_terminals_no_end(self):  # TODO: move to integration or mock muscle
        residues = orderfinder.extract_cterminus(self.data_dir, self.features, "")
        assert residues == {'STAUR_3972': 'ES', 'STAUR_3982': 'GK',
                            'STAUR_3983': 'DS', 'STAUR_3984': 'DS',
                            'STAUR_3985': 'DS'}

    def test_c_terminals_with_end(self):  # TODO: move to integration or mock muscle
        residues = orderfinder.extract_cterminus(self.data_dir, self.features, "STAUR_3982")
        assert residues == {'STAUR_3972': 'ES', 'STAUR_3983': 'DS',
                            'STAUR_3984': 'DS', 'STAUR_3985': 'DS'}


class TestNTerminalExtract(unittest.TestCase):
    def setUp(self):
        # 60 long, only the first 50 should be used
        self.seqs = {"STAUR_3972": "MQPLEGRFAGRTVVVTGAGAGIGHATASRLMREGARVVASDIAQDRLAALEAESPRGALV",
                     "STAUR_3982": "MSQPENEYLSRLRNAVVALREMQQEIDALNHARTEPIAIVGMGCRFPGGASTPEAFWKLL",
                     "STAUR_3985": "MRQAGSPSSPEALQSLVISLVAARTALPVRSIDVREPLSRHGLDSAGAMGLLAELSADLG",
                     "STAUR_3983": "MSVSEADYIARLRKAAITLKEMEGKLGALERARTEPIAIIGMGCRLPGGASTPEAFWKLL",
                     "STAUR_3984": "MNDASSMSTVKRALLAVQEMKARLDAVTRAQTEPIAIIGLGCRLPGGASTPEAFWKLIES"}

        self.data_dir = path.get_full_path(os.path.dirname(__file__), "data", "terminals")

        self.features = [DummyCDS(1, 200, translation=seq, locus_tag=seq_id) for seq_id, seq in self.seqs.items()]
        self.features_by_id = {feature.locus_tag: feature for feature in self.features}

    def test_n_terminals_no_start(self):  # TODO: move to integration or mock muscle
        residues = orderfinder.extract_nterminus(self.data_dir, self.features, "")
        assert residues == {'STAUR_3972': 'L-', 'STAUR_3982': 'ER',
                            'STAUR_3983': 'DK', 'STAUR_3984': 'SQ',
                            'STAUR_3985': 'SV'}

    def test_n_terminals_with_start(self):  # TODO: move to integration or mock muscle
        residues = orderfinder.extract_nterminus(self.data_dir, self.features, "STAUR_3982")
        assert residues == {'STAUR_3972': 'L-', 'STAUR_3983': 'DK',
                            'STAUR_3984': 'SQ', 'STAUR_3985': 'SV'}


class TestOrdering(unittest.TestCase):
    def run_ranking_as_genes(self, n_terms, c_terms, orders):
        genes = {name: DummyCDS(1, 2, locus_tag=name) for name in orders[0]}
        gene_orders = [[genes[k] for k in order] for order in orders]
        res = orderfinder.rank_biosynthetic_orders(n_terms, c_terms, gene_orders)
        if isinstance(orders[0], str):
            return "".join(gene.locus_tag for gene in res)
        return [gene.locus_tag for gene in res]

    def test_permutations_no_start_no_end(self):
        perms = orderfinder.find_possible_orders(list("abc"), "", "")
        assert perms == [["a", "b", "c"], ["a", "c", "b"],
                         ["b", "a", "c"], ["b", "c", "a"],
                         ["c", "a", "b"], ["c", "b", "a"]]

    def test_permutations_end_only(self):
        perms = orderfinder.find_possible_orders(list("abc"), "", "a")
        assert perms == [["b", "c", "a"], ["c", "b", "a"]]

    def test_permutations_start_only(self):
        perms = orderfinder.find_possible_orders(list("abc"), "a", "")
        assert perms == [["a", "b", "c"], ["a", "c", "b"]]

    def test_permutations_start_and_end(self):
        perms = orderfinder.find_possible_orders(list("abc"), "a", "b")
        assert perms == [["a", "c", "b"]]
        perms = orderfinder.find_possible_orders(list("abce"), "a", "b")
        assert perms == [["a", "c", "e", "b"], ["a", "e", "c", "b"]]

    def test_finding_end_gene(self):
        inputs = {"STAUR_3972": ['PKS_KR'],
                  "STAUR_3982": ['PKS_KS', 'PKS_AT', 'PKS_DH2', 'PKS_KR', 'ACP', 'Thioesterase'],
                  "STAUR_3983": ['PKS_KS', 'PKS_AT', 'PKS_DH', 'PKS_KR', 'ACP'],
                  "STAUR_3984": ['PKS_KS', 'PKS_AT', 'PKS_DH', 'PKS_KR', 'ACP']}
        genes = {}
        for name, domains in inputs.items():
            cds = DummyCDS(1, 2, locus_tag=name)
            cds.nrps_pks.domain_names = domains
            genes[name] = cds
        start, end = orderfinder.find_first_and_last_genes(genes.values())
        assert not start
        assert end.get_name() == "STAUR_3982"
        genes["STAUR_3983"].nrps_pks.domain_names.append("TD")
        start, end = orderfinder.find_first_and_last_genes(genes.values())
        assert not start
        assert not end
        genes["STAUR_3984"].nrps_pks.domain_names.append("Thiosterase")
        start, end = orderfinder.find_first_and_last_genes(genes.values())
        assert not start
        assert not end

    def test_finding_start_gene(self):
        inputs = {"STAUR_3972": ['PKS_KR'],
                  "STAUR_3983": ['PKS_KS', 'PKS_AT', 'PKS_DH', 'PKS_KR', 'ACP'],
                  "STAUR_3984": ['PKS_KS', 'PKS_AT', 'PKS_DH', 'PKS_KR', 'ACP'],
                  "STAUR_3985": ['ACP', 'PKS_KS', 'PKS_AT', 'PKS_DH', 'PKS_KR', 'ACP']}
        genes = {}
        for name, domains in inputs.items():
            cds = DummyCDS(1, 2, locus_tag=name)
            cds.nrps_pks.domain_names = domains
            genes[name] = cds
        # no starts
        start, end = orderfinder.find_first_and_last_genes(genes.values())
        assert not start
        assert not end
        # fallback start
        genes["STAUR_3983"].nrps_pks.domain_names = ["PKS_KS", "PKS_AT", "ACP"]
        start, end = orderfinder.find_first_and_last_genes(genes.values())
        assert start.get_name() == "STAUR_3983"
        assert not end
        # two fallback start possibilities
        genes["STAUR_3984"].nrps_pks.domain_names = ["PKS_KS", "PKS_AT", "ACP"]
        start, end = orderfinder.find_first_and_last_genes(genes.values())
        assert not start
        assert not end
        # first-class start
        genes["STAUR_3972"].nrps_pks.domain_names = ["PKS_AT", "ACP"]
        start, end = orderfinder.find_first_and_last_genes(genes.values())
        assert start.get_name() == "STAUR_3972"
        assert not end
        # two possible starts
        genes["STAUR_3984"].nrps_pks.domain_names = ["PKS_AT", "ACP"]
        start, end = orderfinder.find_first_and_last_genes(genes.values())
        assert not start
        assert not end

    def test_pair_polarity(self):
        # first pairing polarity
        n_terms = {"A": "H-", "B": "D-", "C": "H-"}
        c_terms = {"A": "D-", "B": "H-", "C": "H-"}
        possible_orders = ["ABC", "ACB", "BAC", "BCA", "CAB", "CBA"]
        best = self.run_ranking_as_genes(n_terms, c_terms, possible_orders)
        assert best == "ACB"
        # and that the order of orders has no impact
        possible_orders.reverse()
        best = self.run_ranking_as_genes(n_terms, c_terms, possible_orders)
        assert best == "ACB"

        # second pairing polarity
        n_terms = {key: list(reversed(val)) for key, val in n_terms.items()}
        c_terms = {key: list(reversed(val)) for key, val in c_terms.items()}
        possible_orders = ["ABC", "ACB", "BAC", "BCA", "CAB", "CBA"]
        best = self.run_ranking_as_genes(n_terms, c_terms, possible_orders)
        assert best == "ACB"
        # and that the order of orders has no impact
        possible_orders.reverse()
        best = self.run_ranking_as_genes(n_terms, c_terms, possible_orders)
        assert best == "ACB"

    def test_pair_hydrophobic(self):
        # first hydrophobic
        n_terms = {"A": "A-", "B": "X-", "C": "L-"}
        c_terms = {"A": "V-", "B": "I-", "C": "X-"}
        possible_orders = ["ABC", "ACB", "BAC", "BCA", "CAB", "CBA"]
        best = self.run_ranking_as_genes(n_terms, c_terms, possible_orders)
        assert best == "BAC"
        # and that the order of orders has no impact
        possible_orders.reverse()
        best = self.run_ranking_as_genes(n_terms, c_terms, possible_orders)
        assert best == "BAC"

        # second hydrophobic
        n_terms = {key: list(reversed(val)) for key, val in n_terms.items()}
        c_terms = {key: list(reversed(val)) for key, val in c_terms.items()}
        possible_orders = ["ABC", "ACB", "BAC", "BCA", "CAB", "CBA"]
        best = self.run_ranking_as_genes(n_terms, c_terms, possible_orders)
        assert best == "BAC"
        # and that the order of orders has no impact
        possible_orders.reverse()
        best = self.run_ranking_as_genes(n_terms, c_terms, possible_orders)
        assert best == "BAC"

    def test_order_C002271_c19(self):  # pylint: disable=invalid-name
        n_terms = {'STAUR_3972': 'L-', 'STAUR_3982': 'ER',
                   'STAUR_3983': 'DK', 'STAUR_3984': 'SQ',
                   'STAUR_3985': 'SV'}
        c_terms = {'STAUR_3972': 'ES', 'STAUR_3983': 'DS',
                   'STAUR_3984': 'DS', 'STAUR_3985': 'DS'}
        start = ""
        end = "STAUR_3982"
        # there are multiple orders of equal score, but sort for simple testing
        possible_orders = sorted(orderfinder.find_possible_orders(list(n_terms), start, end))
        best = self.run_ranking_as_genes(n_terms, c_terms, possible_orders)
        assert best == ['STAUR_3983', 'STAUR_3972', 'STAUR_3984', 'STAUR_3985', 'STAUR_3982']


class TestEnzymeCounter(unittest.TestCase):
    def run_finder(self, names, all_domains, types=None):
        genes = [DummyCDS(1, 2, locus_tag=name) for name in names]
        for gene in genes:
            gene.nrps_pks.domain_names = all_domains[gene.get_name()]
            if not types:
                gene.nrps_pks.type = domain_identification.classify_feature(all_domains[gene.get_name()])
            else:
                gene.nrps_pks.type = types[gene.get_name()]
        return orderfinder.find_cluster_modular_enzymes(genes)

    def test_C002271_c19(self):  # pylint: disable=invalid-name
        gene_names = ['STAUR_3972', 'STAUR_3982', 'STAUR_3983', 'STAUR_3984', 'STAUR_3985']
        gene_domains = {'STAUR_3985': ['ACP', 'PKS_KS', 'PKS_AT', 'PKS_DH', 'PKS_KR', 'ACP'],
                        'STAUR_3984': ['PKS_KS', 'PKS_AT', 'PKS_DH', 'PKS_KR', 'ACP'],
                        'STAUR_3983': ['PKS_KS', 'PKS_AT', 'PKS_DH', 'PKS_KR', 'ACP'],
                        'STAUR_3982': ['PKS_KS', 'PKS_AT', 'PKS_DH2', 'PKS_KR', 'ACP', 'Thioesterase'],
                        'STAUR_3972': ['PKS_KR']}
        gene_types = {'STAUR_3985': 'PKS-like protein', 'STAUR_3984': 'PKS-like protein',
                      'STAUR_3983': 'PKS-like protein', 'STAUR_3982': 'PKS-like protein',
                      'STAUR_3972': 'other'}
        result = self.run_finder(gene_names, gene_domains, gene_types)
        assert result == (4, 0, 0)

    def test_none(self):
        assert orderfinder.find_cluster_modular_enzymes([]) == (0, 0, 0)

    def test_blended(self):
        names = list("BC")
        domains = {"B": ["PKS_AT"], "C": ["PKS_KS", "PKS_AT"]}
        assert self.run_finder(names, domains) == (2, 0, 0)

        domains = {"B": ["AMP-binding"], "C": ["PKS_AT"]}
        assert self.run_finder(names, domains) == (1, 1, 0)

        domains = {"B": ["AMP-binding", "PKS_AT"], "C": ["PKS_AT"]}
        assert self.run_finder(names, domains) == (1, 0, 0)

        domains = {"B": ["AMP-binding", "PKS_AT", "PKS_KS"], "C": ["PKS_AT"]}
        assert self.run_finder(names, domains) == (1, 0, 0)

        domains = {"B": ["AMP-binding", "Condensation_Dual", "PKS_AT"], "C": ["PKS_AT"]}
        assert self.run_finder(names, domains) == (1, 0, 1)

        domains = {"B": ["AMP-binding"], "C": ["PKS_AT"]}
        assert self.run_finder(names, domains) == (1, 1, 0)
