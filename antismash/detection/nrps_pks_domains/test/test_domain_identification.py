# License: GNU Affero General Public License v3 or later
# A copy of GNU AGPL v3 should have been included in this software package in LICENSE.txt.

# for test files, silence irrelevant and noisy pylint warnings
# pylint: disable=no-self-use,protected-access,missing-docstring

import unittest

from antismash.common.hmmscan_refinement import HMMResult
from antismash.common.test.helpers import DummyRecord, DummyCDS
from antismash.detection.nrps_pks_domains import domain_identification


def dummy_hmm(hit_id="dummy", start=1):
    return HMMResult(hit_id, start, start + 40, 1e-5, 10)


class TestTerminalRemoval(unittest.TestCase):
    def setUp(self):
        self.func = domain_identification.filter_nonterminal_docking_domains
        self.cds = DummyCDS(0, 200)
        self.cds.translation = "A" * 200
        self.record = DummyRecord(features=[self.cds])

    def test_no_terminals(self):
        hits = [dummy_hmm(start=5 + i * 10) for i in range(10)]
        results = self.func(self.record, {self.cds.locus_tag: hits})
        assert results[self.cds.locus_tag] == hits
        # make sure we're not altering in place
        assert id(results[self.cds.locus_tag]) != id(hits)

    def test_single_leading_terminal_in_range(self):
        hits = [dummy_hmm(start=5 + i * 10) for i in range(2)]
        for name in {'NRPS-COM_Nterm', 'NRPS-COM_Cterm', 'PKS_Docking_Cterm',
                     'PKS_Docking_Nterm'}:
            hits[0].hit_id = name
            results = self.func(self.record, {self.cds.locus_tag: hits})
            assert results[self.cds.locus_tag] == hits

    def test_single_terminal_in_mid(self):
        hits = [dummy_hmm(start=50 + i * 10) for i in range(2)]
        for name in {'NRPS-COM_Nterm', 'NRPS-COM_Cterm', 'PKS_Docking_Cterm',
                     'PKS_Docking_Nterm'}:
            hits[0].hit_id = name
            results = self.func(self.record, {self.cds.locus_tag: hits})
            assert results[self.cds.locus_tag] == hits[1:]

    def test_trailing_terminal(self):
        hits = [dummy_hmm(start=5 + i * 50) for i in range(4)]
        print(hits)
        for name in {'NRPS-COM_Nterm', 'NRPS-COM_Cterm', 'PKS_Docking_Cterm',
                     'PKS_Docking_Nterm'}:
            hits[-1].hit_id = name
            results = self.func(self.record, {self.cds.locus_tag: hits})
            print(results)
            assert results[self.cds.locus_tag] == hits


class TestCDSClassification(unittest.TestCase):
    def func(self, domains, ks_subtypes=None):
        return domain_identification.classify_cds(domains, ks_subtypes or [])

    def test_other(self):
        assert self.func(["notPKS", "notNRPS"]) == "other"

    def test_nrps_like(self):
        assert self.func(["AMP-binding", "PKS_KS"]) == "PKS/NRPS-like protein"
        assert self.func(["Condensation_Starter", "PKS_KS"]) == "PKS/NRPS-like protein"

    def test_hybrid(self):
        assert self.func(["AMP-binding", "Condensation_LCL", "PKS_KS"]) == "Hybrid PKS-NRPS"

    def test_nrps(self):
        assert self.func(["AMP-binding", "Condensation_LCL"]) == "NRPS"

    def test_trans_at(self):
        trans = "Type I Trans-AT PKS"
        base = ["PKS_KS", "Trans-AT_docking"]
        assert self.func(base) != trans
        assert self.func(base + ["PKS_KS"], ["Trans-AT-KS", "Modular-KS"]) != trans
        assert self.func(base, ["Trans-AT-KS"]) == trans
        assert self.func(base + ["PKS_AT"], ["Trans-AT-KS"]) != trans

    def test_iterative(self):
        iterative = "Type I Iterative PKS"
        # missing AT
        assert self.func(["PKS_KS"], ["Iterative-KS"]) != iterative
        # too many KS domains
        assert self.func(["PKS_KS"]*3, ["Iterative-KS"]*3) != iterative
        # iterative isn't most common
        assert self.func(["PKS_AT"] + ["PKS_KS"]*2, ["Iterative-KS", "Modular-KS"]) != iterative

        assert self.func(["PKS_KS", "PKS_AT"], ["Iterative-KS"]) == iterative
        assert self.func(["PKS_AT"] + ["PKS_KS"]*2, ["Iterative-KS"]*2) == iterative

    def test_enediyne(self):
        ene = "Type I Enediyne PKS"
        # missing AT
        assert self.func(["PKS_KS"], ["Enediyne-KS"]) != ene
        # too many KS domains
        assert self.func(["PKS_AT"] + ["PKS_KS"]*3, ["Enediyne-KS"]*3) != ene
        # ene isn't most common
        assert self.func(["PKS_AT"] + ["PKS_KS"]*2, ["Enediyne-KS", "Modular-KS"]) != ene

        assert self.func(["PKS_KS", "PKS_AT"], ["Enediyne-KS"]) == ene
        assert self.func(["PKS_AT"] + ["PKS_KS"]*2, ["Enediyne-KS"]*2) == ene

    def test_modular(self):
        mod = "Type I Modular PKS"
        # missing AT
        assert self.func(["PKS_KS"], ["Enediyne-KS"]) != mod
        # more than three KS, regardless of type
        assert self.func(["PKS_AT"] + ["PKS_KS"]*4, ["Enediyne-KS"]*4) == mod
        # three KS, but mod is most
        assert self.func(["PKS_AT"] + ["PKS_KS"]*3, ["Modular-KS"]*3) == mod

    def test_glycopeptide(self):
        glyc = "Glycopeptide NRPS"
        base = ["Cglyc", "Epimerization", "AMP-binding"]
        assert self.func(base) == glyc
        assert self.func(base + ["PKS_KS"]) != glyc
        assert self.func(base + ["PKS_AT"]) != glyc

        # all three domains have to be present
        assert self.func(base[1:]) != glyc
        assert self.func(base[::2]) != glyc
        assert self.func(base[:-1]) != glyc


class TestKSCounter(unittest.TestCase):
    def setUp(self):
        self.counter = domain_identification.KetosynthaseCounter
        self.one_each = ["PKS_KS", "Trans-AT-KS", "Modular-KS", "Enediyne-KS", "Iterative-KS"]

    def test_no_count(self):
        counter = self.counter(["dummy"] * 5)
        for count in vars(counter).values():
            assert count == 0
        assert not counter.modular_is_greatest()
        assert not counter.ene_is_greatest()
        assert not counter.iterative_is_greatest()
        assert not counter.trans_is_greatest()

    def test_single_each(self):
        counter = self.counter(self.one_each)
        for count in vars(counter).values():
            assert count == 1
        assert not counter.modular_is_greatest()
        assert not counter.ene_is_greatest()
        assert not counter.iterative_is_greatest()
        assert not counter.trans_is_greatest()

    def test_specific(self):
        for var_name, ks_name in zip(["trans_at", "modular", "enediyne", "iterative"],
                            ["Trans-AT-KS", "Modular-KS", "Enediyne-KS", "Iterative-KS"]):
            counter = self.counter(self.one_each + [ks_name])
            for key, val in vars(counter).items():
                if key == var_name:
                    assert val == 2
                else:
                    assert val == 1, "%s %s" % (ks_name, vars(counter))
            assert counter.modular_is_greatest() == (ks_name == "Modular-KS")
            assert counter.ene_is_greatest() == (ks_name == "Enediyne-KS")
            assert counter.iterative_is_greatest() == (ks_name == "Iterative-KS")
            assert counter.trans_is_greatest() == (ks_name == "Trans-AT-KS")
