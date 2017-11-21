# License: GNU Affero General Public License v3 or later
# A copy of GNU AGPL v3 should have been included in this software package in LICENSE.txt.

# for test files, silence irrelevant and noisy pylint warnings
# pylint: disable=no-self-use,protected-access,missing-docstring

import unittest

from antismash.common.test import helpers
from antismash.modules.nrps_pks import parsers


class TestNRPSParserMonomerModification(unittest.TestCase):
    def setUp(self):
        self.genes = []
        self.clusters = []
        domain_names = self.gen_domain_names()
        for product in ['not_atpks', 'transatpks']:
            cluster = helpers.DummyCluster(1, 2)
            cluster.products = [product]
            for i in range(7):
                locus_tag = chr(ord('a') + i)
                if i == 6:
                    locus_tag = "all"
                cds = helpers.DummyCDS(1, 2, locus_tag=locus_tag)
                cds.product = product
                cds.nrps_pks.domain_names = domain_names[locus_tag]
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
            res[chr(key)] = [v]
            res[chr(key + len(possible_values))] = [v]*3
        res['all'] = possible_values*3
        return res

    def gen_predictions(self, gene_name, gene_domains, pred):
        # predictions per domain, not gene
        res = {}
        for pks_type in ["KS", "AT"]:
            counter = 0
            for domain in gene_domains:
                if domain == "PKS_" + pks_type:
                    counter += 1
                    res["{}_{}{}".format(gene_name, pks_type, counter)] = pred
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
            ({('all_AT3', 'redmal'), ('all_AT1', 'redmal'), ('all_AT2', 'redmal')},
             {('all_KS2', 'redmal'), ('all_KS3', 'redmal'), ('all_KS1', 'redmal')}),  # ccmal
            ({('all_AT3', 'redemal'), ('all_AT1', 'redemal'), ('all_AT2', 'redemal')},
             {('all_KS2', 'redemal'), ('all_KS3', 'redemal'), ('all_KS1', 'redemal')}),  # ohemal
            ({('all_AT3', 'redmxmal'), ('all_AT1', 'redmxmal'), ('all_AT2', 'redmxmal')},
             {('all_KS2', 'redmxmal'), ('all_KS3', 'redmxmal'), ('all_KS1', 'redmxmal')}),  # ohmxmal
            ({('all_AT3', 'redmmal'), ('all_AT1', 'redmmal'), ('all_AT2', 'redmmal')},
             {('all_KS2', 'redmmal'), ('all_KS3', 'redmmal'), ('all_KS1', 'redmmal')}),  # ohmmal
            ({('all_AT3', 'redmmal'), ('all_AT1', 'redmmal'), ('all_AT2', 'redmmal')},
             {('all_KS2', 'redmmal'), ('all_KS3', 'redmmal'), ('all_KS1', 'redmmal')}),  # ccmmal
            ({('all_AT3', 'redemal'), ('all_AT1', 'redemal'), ('all_AT2', 'redemal')},
             {('all_KS2', 'redemal'), ('all_KS1', 'redemal')}),  # emal
            (set(), set()),  # redmmal
            ({('all_AT3', 'redmmal'), ('all_AT1', 'redmmal'), ('all_AT2', 'redmmal')},
             {('all_KS2', 'redmmal'), ('all_KS1', 'redmmal')}),  # mmal
            ({('all_AT3', 'redmxmal'), ('all_AT1', 'redmxmal'), ('all_AT2', 'redmxmal')},
             {('all_KS2', 'redmxmal'), ('all_KS3', 'redmxmal'), ('all_KS1', 'redmxmal')}),  # ccmxmal
            ({('all_AT3', 'redmxmal'), ('all_AT1', 'redmxmal'), ('all_AT2', 'redmxmal')},
             {('all_KS2', 'redmxmal'), ('all_KS1', 'redmxmal')}),  # mxmal
            (set(), set()),  # redemal
            ({('all_AT3', 'redmal'), ('all_AT1', 'redmal'), ('all_AT2', 'redmal')},
             {('all_KS2', 'redmal'), ('all_KS3', 'redmal'), ('all_KS1', 'redmal')}),  # ohmal
            ({('all_AT3', 'redmal'), ('all_AT1', 'redmal'), ('all_AT2', 'redmal')},
             {('all_KS2', 'redmal'), ('all_KS1', 'redmal')}),   # mal
            ({('all_AT3', 'redemal'), ('all_AT1', 'redemal'), ('all_AT2', 'redemal')},
             {('all_KS2', 'redemal'), ('all_KS3', 'redemal'), ('all_KS1', 'redemal')})]  # ccemal

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
