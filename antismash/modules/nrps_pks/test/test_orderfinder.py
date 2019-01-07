# License: GNU Affero General Public License v3 or later
# A copy of GNU AGPL v3 should have been included in this software package in LICENSE.txt.

# for test files, silence irrelevant and noisy pylint warnings
# pylint: disable=no-self-use,protected-access,missing-docstring

import os
import unittest

from antismash.common import path, secmet
from antismash.common.test.helpers import DummyCDS
from antismash.detection.nrps_pks_domains import domain_identification
from antismash.modules.nrps_pks import orderfinder


class DummyNRPSQualfier(secmet.qualifiers.NRPSPKSQualifier):  # pylint: disable=abstract-method
    def __init__(self):
        super().__init__(strand=1)

    @property
    def domain_names(self):
        return self._domain_names

    @domain_names.setter
    def domain_names(self, names):
        self._domain_names = names


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
        residues = orderfinder.extract_cterminus(self.data_dir, self.features, self.features_by_id["STAUR_3982"])
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
        residues = orderfinder.extract_nterminus(self.data_dir, self.features, None)
        assert residues == {'STAUR_3972': 'L-', 'STAUR_3982': 'ER',
                            'STAUR_3983': 'DK', 'STAUR_3984': 'SQ',
                            'STAUR_3985': 'SV'}

    def test_n_terminals_with_start(self):  # TODO: move to integration or mock muscle
        residues = orderfinder.extract_nterminus(self.data_dir, self.features, self.features_by_id["STAUR_3982"])
        assert residues == {'STAUR_3972': 'L-', 'STAUR_3983': 'DK',
                            'STAUR_3984': 'SQ', 'STAUR_3985': 'SV'}


class TestOrdering(unittest.TestCase):
    def setUp(self):
        self.gene_mapping = {"a": DummyCDS(1, 2, locus_tag="a"),
                             "b": DummyCDS(3, 4, locus_tag="b"),
                             "c": DummyCDS(5, 6, locus_tag="c"),
                             "e": DummyCDS(7, 8, locus_tag="e")}
        self.genes = list(self.gene_mapping.values())

    def run_ordering_simple(self, start, end, gene_list=None):
        start_gene = self.gene_mapping.get(start)
        end_gene = self.gene_mapping.get(end)
        if gene_list is None:
            gene_list = "abc"
        genes = [self.gene_mapping[name] for name in gene_list]
        orders = orderfinder.find_possible_orders(genes, start_gene, end_gene)
        simple_orders = []
        for order in orders:
            simple_orders.append([g.locus_tag for g in order])
        return simple_orders

    def run_ranking_as_genes(self, n_terms, c_terms, orders):
        genes = {name: DummyCDS(locus_tag=name) for name in orders[0]}
        gene_orders = [[genes[k] for k in order] for order in orders]
        res = orderfinder.rank_biosynthetic_orders(n_terms, c_terms, gene_orders)
        return "".join(gene.locus_tag for gene in res)

    def test_permutations_no_start_no_end(self):
        perms = self.run_ordering_simple("", "")
        assert perms == [["a", "b", "c"], ["a", "c", "b"],
                         ["b", "a", "c"], ["b", "c", "a"],
                         ["c", "a", "b"], ["c", "b", "a"]]

    def test_permutations_end_only(self):
        perms = self.run_ordering_simple("", "a")
        assert perms == [["b", "c", "a"], ["c", "b", "a"]]

    def test_permutations_start_only(self):
        perms = self.run_ordering_simple("a", "")
        assert perms == [["a", "b", "c"], ["a", "c", "b"]]

    def test_permutations_start_and_end(self):
        perms = self.run_ordering_simple("a", "b")
        assert perms == [["a", "c", "b"]]
        perms = self.run_ordering_simple("a", "b", "abce")
        assert perms == [["a", "c", "e", "b"], ["a", "e", "c", "b"]]

    def test_same_start_and_end(self):
        with self.assertRaisesRegex(AssertionError, "Using same gene for start and end of ordering"):
            self.run_ordering_simple("a", "a")

    def test_finding_end_gene(self):
        inputs = {"STAUR_3972": ['PKS_KR'],
                  "STAUR_3982": ['PKS_KS', 'PKS_AT', 'PKS_DH2', 'PKS_KR', 'ACP', 'Thioesterase'],
                  "STAUR_3983": ['PKS_KS', 'PKS_AT', 'PKS_DH', 'PKS_KR', 'ACP'],
                  "STAUR_3984": ['PKS_KS', 'PKS_AT', 'PKS_DH', 'PKS_KR', 'ACP']}
        genes = {}
        for name, domains in inputs.items():
            cds = DummyCDS(locus_tag=name)
            cds.nrps_pks = DummyNRPSQualfier()
            cds.nrps_pks.domain_names = domains
            genes[name] = cds
        start, end = orderfinder.find_first_and_last_cds(genes.values())
        assert not start
        assert end.get_name() == "STAUR_3982"
        genes["STAUR_3983"].nrps_pks.domain_names.append("TD")
        start, end = orderfinder.find_first_and_last_cds(genes.values())
        assert not start
        assert not end
        genes["STAUR_3984"].nrps_pks.domain_names.append("Thiosterase")
        start, end = orderfinder.find_first_and_last_cds(genes.values())
        assert not start
        assert not end

    def test_finding_start_gene(self):
        inputs = {"STAUR_3972": ['PKS_KR'],
                  "STAUR_3983": ['PKS_KS', 'PKS_AT', 'PKS_DH', 'PKS_KR', 'ACP'],
                  "STAUR_3984": ['PKS_KS', 'PKS_AT', 'PKS_DH', 'PKS_KR', 'ACP'],
                  "STAUR_3985": ['ACP', 'PKS_KS', 'PKS_AT', 'PKS_DH', 'PKS_KR', 'ACP']}
        genes = {}
        for name, domains in inputs.items():
            cds = DummyCDS(locus_tag=name)
            cds.nrps_pks = DummyNRPSQualfier()
            cds.nrps_pks.domain_names = domains
            genes[name] = cds
        # no starts
        start, end = orderfinder.find_first_and_last_cds(genes.values())
        assert not start
        assert not end
        # fallback start
        genes["STAUR_3983"].nrps_pks.domain_names = ["PKS_KS", "PKS_AT", "ACP"]
        start, end = orderfinder.find_first_and_last_cds(genes.values())
        assert start.get_name() == "STAUR_3983"
        assert not end
        # two fallback start possibilities
        genes["STAUR_3984"].nrps_pks.domain_names = ["PKS_KS", "PKS_AT", "ACP"]
        start, end = orderfinder.find_first_and_last_cds(genes.values())
        assert not start
        assert not end
        # first-class start
        genes["STAUR_3972"].nrps_pks.domain_names = ["PKS_AT", "ACP"]
        start, end = orderfinder.find_first_and_last_cds(genes.values())
        assert start.get_name() == "STAUR_3972"
        assert not end
        # two possible starts
        genes["STAUR_3984"].nrps_pks.domain_names = ["PKS_AT", "ACP"]
        start, end = orderfinder.find_first_and_last_cds(genes.values())
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
        cdss = {}
        for i, name in enumerate(["STAUR_3972", "STAUR_3982", "STAUR_3983",
                                  "STAUR_3984", "STAUR_3985"]):
            cdss[name] = DummyCDS(start=i*10, end=i*10+1, locus_tag=name)

        n_terms = {'STAUR_3972': 'L-', 'STAUR_3982': 'ER',
                   'STAUR_3983': 'DK', 'STAUR_3984': 'SQ',
                   'STAUR_3985': 'SV'}
        c_terms = {'STAUR_3972': 'ES', 'STAUR_3982': '--',
                   'STAUR_3983': 'DS', 'STAUR_3984': 'DS',
                   'STAUR_3985': 'DS'}
        start = None
        end = cdss["STAUR_3982"]
        # there are multiple orders of equal score, but sort for simple testing
        possible_orders = orderfinder.find_possible_orders(list(cdss.values()), start, end)
        best = orderfinder.rank_biosynthetic_orders(n_terms, c_terms, possible_orders)
        best = [gene.get_name() for gene in best]
        assert best == ['STAUR_3983', 'STAUR_3972', 'STAUR_3984', 'STAUR_3985', 'STAUR_3982']

    def test_order_finding_size(self):
        cdss = [DummyCDS() for i in range(11)]
        with self.assertRaisesRegex(AssertionError, "input too large"):
            orderfinder.find_possible_orders(cdss, None, None)


class TestEnzymeCounter(unittest.TestCase):
    def run_finder(self, names, all_domains, types=None):
        genes = [DummyCDS(1, 2, locus_tag=name) for name in names]
        for gene in genes:
            gene.nrps_pks = DummyNRPSQualfier()
            gene.nrps_pks.domain_names = all_domains[gene.get_name()]
            if not types:
                gene.nrps_pks.type = domain_identification.classify_cds(all_domains[gene.get_name()])
            else:
                gene.nrps_pks.type = types[gene.get_name()]
        results = orderfinder.find_supercluster_modular_enzymes(genes)
        return ([cds.get_name() for cds in results[0]], results[1], results[2])

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
        expected_pks = ["STAUR_3982", "STAUR_3983", "STAUR_3984", "STAUR_3985"]
        assert result == (expected_pks, 0, 0)

    def test_none(self):
        assert orderfinder.find_supercluster_modular_enzymes([]) == ([], 0, 0)

    def test_blended(self):
        names = list("BC")
        domains = {"B": ["PKS_AT"], "C": ["PKS_KS", "PKS_AT"]}
        assert self.run_finder(names, domains) == (["B", "C"], 0, 0)

        domains = {"B": ["AMP-binding"], "C": ["PKS_AT"]}
        assert self.run_finder(names, domains) == (["C"], 1, 0)

        domains = {"B": ["AMP-binding", "PKS_AT"], "C": ["PKS_AT"]}
        assert self.run_finder(names, domains) == (["C"], 0, 0)

        domains = {"B": ["AMP-binding", "PKS_AT", "PKS_KS"], "C": ["PKS_AT"]}
        assert self.run_finder(names, domains) == (["C"], 0, 0)

        domains = {"B": ["AMP-binding", "Condensation_Dual", "PKS_AT"], "C": ["PKS_AT"]}
        assert self.run_finder(names, domains) == (["C"], 0, 1)

        domains = {"B": ["AMP-binding"], "C": ["PKS_AT"]}
        assert self.run_finder(names, domains) == (["C"], 1, 0)
