# License: GNU Affero General Public License v3 or later
# A copy of GNU AGPL v3 should have been included in this software package in LICENSE.txt.

# for test files, silence irrelevant and noisy pylint warnings
# pylint: disable=no-self-use,protected-access,missing-docstring

import unittest

from antismash.common.test import helpers
from antismash.modules.nrps_pks import parsers

from .test_orderfinder import DummyNRPSQualfier

class TestNRPSParserMonomerModification(unittest.TestCase):
    def setUp(self):
        self.genes = []
        self.clusters = []
        domain_names = self.gen_domain_names()
        for product in ['not_atpks', 'transatpks']:
            cluster = helpers.DummyCluster(1, 2, products=[product])
            assert cluster.products == (product,)
            for i in range(7):
                locus_tag = chr(ord('a') + i)
                if i == 6:
                    locus_tag = "all"
                cds = helpers.DummyCDS(1, 2, locus_tag=locus_tag)
                cds.product = product
                cds.nrps_pks = DummyNRPSQualfier()
                cds.nrps_pks.domain_names = domain_names["nrpspksdomains_" + locus_tag]
                cds.cluster = cluster
                cluster.add_cds(cds)
                self.genes.append(cds)
            self.clusters.append(cluster)
        self.predictions = ['redmxmal', 'ccmal', 'ohemal', 'ohmxmal', 'ohmmal',
                            'ccmmal', 'emal', 'redmmal', 'mmal', 'ccmxmal',
                            'mxmal', 'redemal', 'ohmal', 'mal', 'ccemal']

    def gen_domain_names(self):
        # generates a dict with values containing
        # - each type in a list once alone
        # - each type in a list 3 times with no other type
        # - all types in a list, each type repeated 3 times
        res = {}
        possible_values = ['nomatch', 'PKS_AT', 'PKS_KR', 'PKS_KS', 'PKS_DH', 'PKS_ER']
        for n, v in enumerate(possible_values):
            key = ord('a') + n
            res['nrpspksdomains_' + chr(key)] = [v]
            res['nrpspksdomains_' + chr(key + len(possible_values))] = [v]*3
        res['nrpspksdomains_all'] = possible_values*3
        return res

    def gen_predictions(self, gene_name, gene_domains, pred):
        # predictions per domain, not gene
        res = {}
        for pks_type in ["KS", "AT"]:
            counter = 0
            for domain in gene_domains:
                if domain == "PKS_" + pks_type:
                    counter += 1
                    res["nrpspksdomains_{}_{}{}".format(gene_name, pks_type, counter)] = pred
        return res

    def detect_change(self, preds, original, expected):
        changed = set()
        for domain, new in preds.items():
            if new != original:
                changed.add((domain, new))
        assert changed == expected

    def test_insert_modified_monomers(self):
        # in pairs of (at, trans-at)
        expected_changes = [(set(), set()),  # redmxmal
            ({('nrpspksdomains_all_AT3', 'redmal'), ('nrpspksdomains_all_AT1', 'redmal'), ('nrpspksdomains_all_AT2', 'redmal')},
             {('nrpspksdomains_all_KS2', 'redmal'), ('nrpspksdomains_all_KS3', 'redmal'), ('nrpspksdomains_all_KS1', 'redmal')}),  # ccmal
            ({('nrpspksdomains_all_AT3', 'redemal'), ('nrpspksdomains_all_AT1', 'redemal'), ('nrpspksdomains_all_AT2', 'redemal')},
             {('nrpspksdomains_all_KS2', 'redemal'), ('nrpspksdomains_all_KS3', 'redemal'), ('nrpspksdomains_all_KS1', 'redemal')}),  # ohemal
            ({('nrpspksdomains_all_AT3', 'redmxmal'), ('nrpspksdomains_all_AT1', 'redmxmal'), ('nrpspksdomains_all_AT2', 'redmxmal')},
             {('nrpspksdomains_all_KS2', 'redmxmal'), ('nrpspksdomains_all_KS3', 'redmxmal'), ('nrpspksdomains_all_KS1', 'redmxmal')}),  # ohmxmal
            ({('nrpspksdomains_all_AT3', 'redmmal'), ('nrpspksdomains_all_AT1', 'redmmal'), ('nrpspksdomains_all_AT2', 'redmmal')},
             {('nrpspksdomains_all_KS2', 'redmmal'), ('nrpspksdomains_all_KS3', 'redmmal'), ('nrpspksdomains_all_KS1', 'redmmal')}),  # ohmmal
            ({('nrpspksdomains_all_AT3', 'redmmal'), ('nrpspksdomains_all_AT1', 'redmmal'), ('nrpspksdomains_all_AT2', 'redmmal')},
             {('nrpspksdomains_all_KS2', 'redmmal'), ('nrpspksdomains_all_KS3', 'redmmal'), ('nrpspksdomains_all_KS1', 'redmmal')}),  # ccmmal
            ({('nrpspksdomains_all_AT3', 'redemal'), ('nrpspksdomains_all_AT1', 'redemal'), ('nrpspksdomains_all_AT2', 'redemal')},
             {('nrpspksdomains_all_KS2', 'redemal'), ('nrpspksdomains_all_KS1', 'redemal')}),  # emal
            (set(), set()),  # redmmal
            ({('nrpspksdomains_all_AT3', 'redmmal'), ('nrpspksdomains_all_AT1', 'redmmal'), ('nrpspksdomains_all_AT2', 'redmmal')},
             {('nrpspksdomains_all_KS2', 'redmmal'), ('nrpspksdomains_all_KS1', 'redmmal')}),  # mmal
            ({('nrpspksdomains_all_AT3', 'redmxmal'), ('nrpspksdomains_all_AT1', 'redmxmal'), ('nrpspksdomains_all_AT2', 'redmxmal')},
             {('nrpspksdomains_all_KS2', 'redmxmal'), ('nrpspksdomains_all_KS3', 'redmxmal'), ('nrpspksdomains_all_KS1', 'redmxmal')}),  # ccmxmal
            ({('nrpspksdomains_all_AT3', 'redmxmal'), ('nrpspksdomains_all_AT1', 'redmxmal'), ('nrpspksdomains_all_AT2', 'redmxmal')},
             {('nrpspksdomains_all_KS2', 'redmxmal'), ('nrpspksdomains_all_KS1', 'redmxmal')}),  # mxmal
            (set(), set()),  # redemal
            ({('nrpspksdomains_all_AT3', 'redmal'), ('nrpspksdomains_all_AT1', 'redmal'), ('nrpspksdomains_all_AT2', 'redmal')},
             {('nrpspksdomains_all_KS2', 'redmal'), ('nrpspksdomains_all_KS3', 'redmal'), ('nrpspksdomains_all_KS1', 'redmal')}),  # ohmal
            ({('nrpspksdomains_all_AT3', 'redmal'), ('nrpspksdomains_all_AT1', 'redmal'), ('nrpspksdomains_all_AT2', 'redmal')},
             {('nrpspksdomains_all_KS2', 'redmal'), ('nrpspksdomains_all_KS1', 'redmal')}),   # mal
            ({('nrpspksdomains_all_AT3', 'redemal'), ('nrpspksdomains_all_AT1', 'redemal'), ('nrpspksdomains_all_AT2', 'redemal')},
             {('nrpspksdomains_all_KS2', 'redemal'), ('nrpspksdomains_all_KS3', 'redemal'), ('nrpspksdomains_all_KS1', 'redemal')})]  # ccemal

        for pred, expected in zip(self.predictions, expected_changes):
            for index in [0, 1]:  # AT and transAT
                for gene in self.clusters[index].cds_children:
                    predictions = self.gen_predictions(gene.get_name(), gene.nrps_pks.domain_names, pred)
                    modified = dict(predictions)
                    parsers.modify_monomer_predictions([gene], modified)
                    if gene.get_name() == "all":
                        self.detect_change(modified, pred, expected[index])
                    else:
                        self.detect_change(modified, pred, set())
