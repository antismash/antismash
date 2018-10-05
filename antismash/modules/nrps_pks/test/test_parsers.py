# License: GNU Affero General Public License v3 or later
# A copy of GNU AGPL v3 should have been included in this software package in LICENSE.txt.

# for test files, silence irrelevant and noisy pylint warnings
# pylint: disable=no-self-use,protected-access,missing-docstring

import functools
import unittest

from antismash.common.secmet import Region
from antismash.common.test import helpers
from antismash.modules.nrps_pks import parsers

from .test_orderfinder import DummyNRPSQualfier


class TestNRPSParserMonomerModification(unittest.TestCase):
    def setUp(self):
        self.genes = []
        self.regions = []
        domain_names = self.gen_domain_names()
        for product in ['not_atpks', 'transatpks']:
            cluster = helpers.DummyCluster(1, 2, product=product)
            supercluster = helpers.DummySuperCluster([cluster])
            self.regions.append(Region(superclusters=[supercluster]))
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
                self.regions[-1].add_cds(cds)
                assert cds.region == self.regions[-1]
        self.predictions = ['redmxmal', 'ccmal', 'ohemal', 'ohmxmal', 'ohmmal',
                            'ccmmal', 'emal', 'redmmal', 'mmal', 'ccmxmal',
                            'mxmal', 'redemal', 'ohmal', 'mal', 'ccemal']

    @staticmethod
    def gen_specific_domain_name(locus, pks_variant, counter):
        return "nrpspksdomains_{}_PKS_{}.{}".format(locus, pks_variant, counter)

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
                    res[self.gen_specific_domain_name(gene_name, pks_type, counter)] = pred
        return res

    def detect_change(self, preds, original, expected):
        changed = set()
        for domain, new in preds.items():
            if new != original:
                changed.add((domain, new))
        assert changed == expected

    def test_insert_modified_monomers(self):
        dom_name = functools.partial(self.gen_specific_domain_name, "all")

        # in pairs of (at, trans-at)
        expected_changes = [
            # redmxmal
            (set(), set()),
            # ccmal
            ({(dom_name("AT", 3), 'redmal'), (dom_name("AT", 1), 'redmal'), (dom_name("AT", 2), 'redmal')},
             {(dom_name("KS", 2), 'redmal'), (dom_name("KS", 3), 'redmal'), (dom_name("KS", 1), 'redmal')}),
            # ohemal
            ({(dom_name("AT", 3), 'redemal'), (dom_name("AT", 1), 'redemal'), (dom_name("AT", 2), 'redemal')},
             {(dom_name("KS", 2), 'redemal'), (dom_name("KS", 3), 'redemal'), (dom_name("KS", 1), 'redemal')}),
            # ohmxmal
            ({(dom_name("AT", 3), 'redmxmal'), (dom_name("AT", 1), 'redmxmal'), (dom_name("AT", 2), 'redmxmal')},
             {(dom_name("KS", 2), 'redmxmal'), (dom_name("KS", 3), 'redmxmal'), (dom_name("KS", 1), 'redmxmal')}),
            # ohmmal
            ({(dom_name("AT", 3), 'redmmal'), (dom_name("AT", 1), 'redmmal'), (dom_name("AT", 2), 'redmmal')},
             {(dom_name("KS", 2), 'redmmal'), (dom_name("KS", 3), 'redmmal'), (dom_name("KS", 1), 'redmmal')}),
            # ccmal
            ({(dom_name("AT", 3), 'redmmal'), (dom_name("AT", 1), 'redmmal'), (dom_name("AT", 2), 'redmmal')},
             {(dom_name("KS", 2), 'redmmal'), (dom_name("KS", 3), 'redmmal'), (dom_name("KS", 1), 'redmmal')}),
            # emal
            ({(dom_name("AT", 3), 'redemal'), (dom_name("AT", 1), 'redemal'), (dom_name("AT", 2), 'redemal')},
             {(dom_name("KS", 2), 'redemal'), (dom_name("KS", 1), 'redemal')}),
            # redmmal
            (set(), set()),
            # mmal
            ({(dom_name("AT", 3), 'redmmal'), (dom_name("AT", 1), 'redmmal'), (dom_name("AT", 2), 'redmmal')},
             {(dom_name("KS", 2), 'redmmal'), (dom_name("KS", 1), 'redmmal')}),
            # ccmxmal
            ({(dom_name("AT", 3), 'redmxmal'), (dom_name("AT", 1), 'redmxmal'), (dom_name("AT", 2), 'redmxmal')},
             {(dom_name("KS", 2), 'redmxmal'), (dom_name("KS", 3), 'redmxmal'), (dom_name("KS", 1), 'redmxmal')}),
            # mxmal
            ({(dom_name("AT", 3), 'redmxmal'), (dom_name("AT", 1), 'redmxmal'), (dom_name("AT", 2), 'redmxmal')},
             {(dom_name("KS", 2), 'redmxmal'), (dom_name("KS", 1), 'redmxmal')}),
            # redemal
            (set(), set()),
            # ohmal
            ({(dom_name("AT", 3), 'redmal'), (dom_name("AT", 1), 'redmal'), (dom_name("AT", 2), 'redmal')},
             {(dom_name("KS", 2), 'redmal'), (dom_name("KS", 3), 'redmal'), (dom_name("KS", 1), 'redmal')}),
            # mal
            ({(dom_name("AT", 3), 'redmal'), (dom_name("AT", 1), 'redmal'), (dom_name("AT", 2), 'redmal')},
             {(dom_name("KS", 2), 'redmal'), (dom_name("KS", 1), 'redmal')}),
            # ccemal
            ({(dom_name("AT", 3), 'redemal'), (dom_name("AT", 1), 'redemal'), (dom_name("AT", 2), 'redemal')},
             {(dom_name("KS", 2), 'redemal'), (dom_name("KS", 3), 'redemal'), (dom_name("KS", 1), 'redemal')}),
        ]

        for pred, expected in zip(self.predictions, expected_changes):
            for index in [0, 1]:  # AT and transAT
                for gene in self.regions[index].cds_children:
                    predictions = self.gen_predictions(gene.get_name(), gene.nrps_pks.domain_names, pred)
                    modified = dict(predictions)
                    parsers.modify_monomer_predictions([gene], modified)
                    if gene.get_name() == "all":
                        self.detect_change(modified, pred, expected[index])
                    else:
                        self.detect_change(modified, pred, set())
