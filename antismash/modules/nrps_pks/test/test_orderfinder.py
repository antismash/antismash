# License: GNU Affero General Public License v3 or later
# A copy of GNU AGPL v3 should have been included in this software package in LICENSE.txt.

# for test files, silence irrelevant and noisy pylint warnings
# pylint: disable=use-implicit-booleaness-not-comparison,protected-access,missing-docstring

import unittest

from antismash.common import secmet
from antismash.common.test.helpers import DummyAntismashDomain, DummyCDS
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


class DummyModule(secmet.features.Module):
    def __init__(self, **kwargs):
        domains = [DummyAntismashDomain(domain=dom) for dom in kwargs.pop("domains")]
        super().__init__(domains, module_type=secmet.features.module.ModuleType.PKS, **kwargs)


class TestOrdering(unittest.TestCase):
    def setUp(self):
        self.gene_mapping = {"a": DummyCDS(1, 2, locus_tag="a"),
                             "b": DummyCDS(3, 4, locus_tag="b"),
                             "c": DummyCDS(5, 6, locus_tag="c"),
                             "e": DummyCDS(7, 8, locus_tag="e")}
        self.genes = list(self.gene_mapping.values())

    def run_ordering_simple(self, start, end, gene_list=None, chains=None):
        start_gene = self.gene_mapping.get(start)
        end_gene = self.gene_mapping.get(end)
        if gene_list is None:
            gene_list = "abc"
        genes = [self.gene_mapping[name] for name in gene_list]
        orders = orderfinder.find_possible_orders(genes, start_gene, end_gene, chains or {})
        simple_orders = []
        for order in orders:
            simple_orders.append([g.locus_tag for g in order])
        return simple_orders

    def run_ranking_as_genes(self, n_terms, c_terms, orders):
        genes = {name: DummyCDS(locus_tag=name) for name in orders[0]}
        gene_orders = [[genes[k] for k in order] for order in orders]
        res = orderfinder.rank_biosynthetic_orders(n_terms, c_terms, gene_orders, {})
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
        inputs = {
            "STAUR_3972": {
                "domains": ['PKS_KR'],
                "complete": False
            },
            "STAUR_3982": {
                "domains": ['PKS_KS', 'PKS_AT', 'PKS_DH2', 'PKS_KR', 'ACP', 'Thioesterase'],
                "complete": True,
                "final": True,
            },
            "STAUR_3983": {
                "domains": ['PKS_KS', 'PKS_AT', 'PKS_DH', 'PKS_KR', 'ACP'],
                "complete": True,
            },
        }
        genes = {}
        for name, details in inputs.items():
            cds = DummyCDS(locus_tag=name)
            cds.add_module(DummyModule(**details))
            cds.nrps_pks = DummyNRPSQualfier()
            cds.nrps_pks.domain_names = details["domains"]
            genes[name] = cds
        start, end = orderfinder.find_first_and_last_cds(genes.values())
        assert not start
        assert end.get_name() == "STAUR_3982"

        # set up a duplicate
        assert not genes["STAUR_3983"].modules[-1].is_final_module()
        genes["STAUR_3983"].modules[-1]._is_final = True
        assert genes["STAUR_3983"].modules[-1].is_final_module()
        # make sure the duplicate causes no single CDS to be returned
        start, end = orderfinder.find_first_and_last_cds(genes.values())
        assert not start
        assert not end

    def test_cross_cds_with_terminal(self):
        other = DummyCDS(locus_tag="other")
        normal = DummyModule(domains=["PKS_KS", "PKS_AT", "ACP"], complete=True, final=False)
        other.add_module(normal)

        head = DummyCDS(locus_tag="head")
        tail = DummyCDS(locus_tag="tail")
        split = DummyModule(domains=["PKS_KS", "PKS_AT", "ACP", "Thioesterase"], complete=True, final=True)
        for cds in [head, tail]:
            cds.add_module(split)
        assert split.is_final_module()
        split._parent_cds_names = ["head", "tail"]
        start, end = orderfinder.find_first_and_last_cds([head, tail, other])
        assert not start
        assert end is tail

    def test_finding_start_gene(self):
        inputs = {
            "STAUR_3972": {
                "domains": ['PKS_KR'],
                "complete": False
            },
            "STAUR_3983": {
                "domains": ['PKS_KS', 'PKS_AT', 'PKS_DH', 'PKS_KR', 'ACP'],
                "complete": True,
            },
            "STAUR_3984": {
                "domains": ['PKS_KS', 'PKS_AT', 'PKS_DH', 'PKS_KR', 'ACP'],
                "complete": True,
            },
        }
        genes = {}
        for name, details in inputs.items():
            cds = DummyCDS(locus_tag=name)
            cds.add_module(DummyModule(**details))
            cds.nrps_pks = DummyNRPSQualfier()
            cds.nrps_pks.domain_names = details["domains"]
            genes[name] = cds

        # no starts
        start, end = orderfinder.find_first_and_last_cds(genes.values())
        assert not start
        assert not end

        # fallback start
        for carrier in ["ACP", "PKS_PP"]:
            genes["STAUR_3983"].nrps_pks.domain_names = ["PKS_KS", "PKS_AT", carrier]
            start, end = orderfinder.find_first_and_last_cds(genes.values())
            assert start.get_name() == "STAUR_3983"
            assert not end

        # two fallback start possibilities
        genes["STAUR_3984"].nrps_pks.domain_names = ["PKS_KS", "PKS_AT", "ACP"]
        start, end = orderfinder.find_first_and_last_cds(genes.values())
        assert not start
        assert not end

        # first-class start
        genes["STAUR_3972"].modules[-1]._is_starter = True
        assert genes["STAUR_3972"].modules[-1].is_starter_module()
        start, end = orderfinder.find_first_and_last_cds(genes.values())
        assert start.get_name() == "STAUR_3972"
        assert not end

        # two possible starts
        genes["STAUR_3984"].modules[-1]._is_starter = True
        assert genes["STAUR_3984"].modules[-1].is_starter_module()
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

    def test_order_CP002271_c19(self):  # pylint: disable=invalid-name
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
        possible_orders = orderfinder.find_possible_orders(list(cdss.values()), start, end, {})
        best = orderfinder.rank_biosynthetic_orders(n_terms, c_terms, possible_orders, {})
        best = [gene.get_name() for gene in best]
        assert best == ['STAUR_3983', 'STAUR_3972', 'STAUR_3984', 'STAUR_3985', 'STAUR_3982']

        # and try again with neighbours set up
        neighbours = {
            "STAUR_3983": "STAUR_3982",
            "STAUR_3984": "STAUR_3983",
            "STAUR_3985": "STAUR_3984",
        }
        best = orderfinder.rank_biosynthetic_orders(n_terms, c_terms, possible_orders, neighbours)
        best = [gene.get_name() for gene in best]
        assert best == ['STAUR_3972', 'STAUR_3985', 'STAUR_3984', 'STAUR_3983', 'STAUR_3982']

    def test_order_finding_size(self):
        cdss = [DummyCDS() for i in range(11)]
        with self.assertRaisesRegex(AssertionError, "input too large"):
            orderfinder.find_possible_orders(cdss, None, None, {})

    def test_split_module(self):
        features = {}
        for i, name in enumerate("ABCDEF"):
            features[name] = (DummyCDS(start=i * 10, end=i * 10 + 6, locus_tag=name))
        # ensure empty chains cause no issues
        singles = orderfinder.find_possible_orders(list(features.values()), None, features["E"], {})
        assert len(singles) == 120

        def check(order, names):
            indices = [order.index(features[name]) for name in names]
            # +1 because the slice changes the maths
            return all(index - indices[0] - (i + 1) == 0 for i, index in enumerate(indices[1:]))

        # check chains are respected
        chains = {
            features["A"]: features["B"],
            features["E"]: features["F"],  # this should still be respected, even with E as the 'end'
        }
        chained = orderfinder.find_possible_orders(list(features.values()), None, features["E"], chains)
        assert len(chained) == 6
        # regardless of where A and B land, they should be consecutive and in that order
        assert all(check(order, "AB") for order in chained)
        # the tail of every order should be EF, since E was marked as the 'end' CDS
        assert all(order[-2:] == [features["E"], features["F"]] for order in chained), chained

        # check complex chains
        chains = {
            features["B"]: features["C"],
            features["C"]: features["D"],
            features["D"]: features["E"],
        }
        chained = orderfinder.find_possible_orders(list(features.values()), None, None, chains)
        assert len(chained) == 6
        # again, regardless of where the block is, the order within the block must be fixed
        assert all(check(order, "BCDE") for order in chained)

    def test_chained_end(self):
        # catches an edge case where order finding generated duplicates of genes
        # specifically where a single gene contained a termination domain and was
        # provided as the "end" cds, where it was also the tail end of a chain of
        # cross-CDS modules
        cdses = [DummyCDS(1, 2, locus_tag="A"), DummyCDS(3, 4, locus_tag="B")]
        orders = orderfinder.find_possible_orders(cdses, start_cds=None, end_cds=cdses[-1],
                                                  chains={cdses[0]: cdses[1]})
        assert len(orders) == 1
        order = orders[0]
        assert len(order) == len(set(order))
        assert order == cdses


class TestEnzymeCounter(unittest.TestCase):
    def run_finder(self, names, modules_by_cds):
        genes = [DummyCDS(1, 2, locus_tag=name) for name in names]
        for gene in genes:
            for domains in modules_by_cds[gene.get_name()]:
                gene.add_module(DummyModule(domains=domains, complete=True))
        results = orderfinder.find_candidate_cluster_modular_enzymes(genes)
        return ([cds.get_name() for cds in results[0]], results[1], results[2])

    def test_CP002271_c19(self):  # pylint: disable=invalid-name
        gene_names = ['STAUR_3972', 'STAUR_3982', 'STAUR_3983', 'STAUR_3984', 'STAUR_3985']
        gene_domains = {'STAUR_3985': [('ACP',), ('PKS_KS', 'PKS_AT', 'PKS_DH', 'PKS_KR', 'ACP')],
                        'STAUR_3984': [('PKS_KS', 'PKS_AT', 'PKS_DH', 'PKS_KR', 'ACP')],
                        'STAUR_3983': [('PKS_KS', 'PKS_AT', 'PKS_DH', 'PKS_KR', 'ACP')],
                        'STAUR_3982': [('PKS_KS', 'PKS_AT', 'PKS_DH2', 'PKS_KR', 'ACP', 'Thioesterase')],
                        'STAUR_3972': [('PKS_KR',)]}
        result = self.run_finder(gene_names, gene_domains)
        expected_pks = ["STAUR_3982", "STAUR_3983", "STAUR_3984", "STAUR_3985"]
        assert result == (expected_pks, 0, 0)

    def test_none(self):
        assert orderfinder.find_candidate_cluster_modular_enzymes([]) == ([], 0, 0)

    def test_blended(self):
        names = list("BC")
        pkses = [[(f"PKS_{kind}",)] for kind in ["KS", "AT"]]
        nrpses = [[(kind,)] for kind in ["AMP-binding", "A-OX", "Condensation"]]
        other = [("T",)]

        for pks in pkses:
            modules = {"B": pks, "C": other}
            assert self.run_finder(names, modules) == (["B"], 0, 0)

        for nrps in nrpses:
            modules = {"B": nrps, "C": pkses[-1]}
            assert self.run_finder(names, modules) == (["C"], 1, 0)
            modules = {"B": nrps, "C": other}
            assert self.run_finder(names, modules) == ([], 1, 0)

            for pks in pkses:
                modules = {"B": pks + nrps, "C": pks}
                assert self.run_finder(names, modules) == (["C"], 0, 1)

                modules = {"B": pks + nrps, "C": nrps}
                assert self.run_finder(names, modules) == ([], 1, 1)

                modules = {"B": nrps + pks + other, "C": pks}
                assert self.run_finder(names, modules) == (["C"], 0, 1)


class TestFollowers(unittest.TestCase):
    def test_gaps(self):
        cdses = []
        for i, name in enumerate("ABC"):
            cdses.append(DummyCDS(start=i*10, end=i*10+1, strand=-1 if name == "B" else 1, locus_tag=name))
        assert orderfinder.get_follower_genes(cdses) == {}

    def test_no_gaps_forward(self):
        cdses = []
        for i, name in enumerate("ABC"):
            cdses.append(DummyCDS(start=i*10, end=i*10+1, strand=1, locus_tag=name))
        assert orderfinder.get_follower_genes(cdses) == {
            "A": "B",
            "B": "C",
        }

    def test_no_gaps_reverse(self):
        cdses = []
        for i, name in enumerate("ABC"):
            cdses.append(DummyCDS(start=i*10, end=i*10+1, strand=-1, locus_tag=name))
        assert orderfinder.get_follower_genes(cdses) == {
            "C": "B",
            "B": "A",
        }

    def test_mix(self):
        cdses = [
            DummyCDS(strand=1, locus_tag="A"),
            DummyCDS(strand=1, locus_tag="B"),
            DummyCDS(strand=-1, locus_tag="C"),
            DummyCDS(strand=-1, locus_tag="D"),
            DummyCDS(strand=1, locus_tag="E"),
        ]
        assert orderfinder.get_follower_genes(cdses) == {
            "A": "B",
            "D": "C",
        }
