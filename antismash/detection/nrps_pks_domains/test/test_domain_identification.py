# License: GNU Affero General Public License v3 or later
# A copy of GNU AGPL v3 should have been included in this software package in LICENSE.txt.

# for test files, silence irrelevant and noisy pylint warnings
# pylint: disable=use-implicit-booleaness-not-comparison,protected-access,missing-docstring

import unittest
from unittest.mock import patch

from antismash.common import fasta, subprocessing, utils
from antismash.common.secmet import Record
from antismash.common.test.helpers import DummyCDS, DummyHMMResult, DummyRecord
from antismash.detection.nrps_pks_domains import domain_identification


def dummy_hmm(hit_id="dummy", start=1):
    return DummyHMMResult(hit_id, start, start + 40, 1e-5, 10)


class TestTerminalRemoval(unittest.TestCase):
    def setUp(self):
        self.func = domain_identification.filter_nonterminal_docking_domains
        self.cds = DummyCDS(0, 200)
        self.cds.translation = "A" * 200
        self.record = DummyRecord(features=[self.cds])
        self.terminals = [
            'NRPS-COM_Nterm',
            'NRPS-COM_Cterm',
            'PKS_Docking_Cterm',
            'PKS_Docking_Nterm'
        ]

    def test_no_terminals(self):
        hits = [dummy_hmm(start=5 + i * 10) for i in range(10)]
        results = self.func(self.record, {self.cds.locus_tag: hits})
        assert results[self.cds.locus_tag] == hits
        # make sure we're not altering in place
        assert id(results[self.cds.locus_tag]) != id(hits)

    def test_single_leading_terminal_in_range(self):
        hits = [dummy_hmm(start=5 + i * 10) for i in range(2)]
        for name in self.terminals:
            hits[0]._hit_id = name
            results = self.func(self.record, {self.cds.locus_tag: hits})
            assert results[self.cds.locus_tag] == hits

    def test_single_terminal_in_mid(self):
        hits = [dummy_hmm(start=50 + i * 10) for i in range(2)]
        for name in self.terminals:
            hits[0]._hit_id = name
            results = self.func(self.record, {self.cds.locus_tag: hits})
            assert results[self.cds.locus_tag] == hits[1:]

    def test_trailing_terminal(self):
        hits = [dummy_hmm(start=5 + i * 50) for i in range(4)]
        for name in self.terminals:
            hits[-1]._hit_id = name
            results = self.func(self.record, {self.cds.locus_tag: hits})
            assert results[self.cds.locus_tag] == hits


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


@patch.object(utils, "get_hmm_lengths", return_value={})
@patch.object(subprocessing, "run_hmmscan", return_value={})
@patch.object(domain_identification, "refine_hmmscan_results")
class TestSubtypeFinding(unittest.TestCase):
    def setUp(self):
        self.record = DummyRecord(seq="A" * 100)
        self.record.add_cds_feature(DummyCDS(locus_tag="left", start=0, end=30))
        self.record.add_cds_feature(DummyCDS(locus_tag="right", start=50, end=80))
        self.existing_hits = {
            "left": [
                DummyHMMResult("a", start=0, end=4),
                DummyHMMResult("b", start=5, end=8),
            ],
            "right": [
                DummyHMMResult("b", start=0, end=4),
                DummyHMMResult("c", start=5, end=8),
            ]
        }
    def find(self, target, existing=None, callback=None):
        if existing is None:
            existing = self.existing_hits
        return domain_identification.find_subtypes(target, "dummy/path", existing,
                                                   self.record, modifier_callback=callback)

    def test_no_relevant_hits(self, refine, hmmscan, lengths):
        result = self.find("other")
        assert not result
        assert not refine.called
        assert not hmmscan.called
        assert not lengths.called
        for hits in self.existing_hits.values():
            assert not any(hit.internal_hits for hit in hits)

    def test_locations(self, refine, hmmscan, _lengths):
        new_hits = {
            "left": [
                DummyHMMResult("b_subtype1", start=0, end=3),
                DummyHMMResult("b_subtype2", start=6, end=7),
            ],
            "right": [
                DummyHMMResult("b_subtype1", start=6, end=8),
            ],
        }
        # set the mock return values
        refine.return_value = new_hits
        # exact value doesn't matter, as long as it's not empty
        hmmscan.return_value = new_hits["left"]

        results = self.find("b")
        # only the second hit matches in "left"
        expected_hit = new_hits["left"][1]
        assert results == {"left": [expected_hit]}
        assert not self.existing_hits["left"][0].internal_hits
        assert self.existing_hits["left"][1].internal_hits == (expected_hit,)

    def test_multiple_per_gene(self, refine, hmmscan, _lengths):
        self.existing_hits["left"][1] = DummyHMMResult("a", start=5, end=8)
        new_hits = {
            "left": [
                DummyHMMResult("a_subtype1", start=0, end=3),
                DummyHMMResult("a_subtype1", start=6, end=7),
            ],
        }
        # set the mock return values
        refine.return_value = new_hits
        # exact value doesn't matter, as long as it's not empty
        hmmscan.return_value = new_hits["left"]

        results = self.find("a")
        # only the second hit matches in "left"
        assert results == new_hits
        for original, new_hit in zip(self.existing_hits["left"], new_hits["left"]):
            assert original.internal_hits == (new_hit,)


@patch.object(fasta, "get_fasta_from_features", return_value={})
@patch.object(domain_identification, "find_ab_motifs", return_value={})
class TestModuleMerging(unittest.TestCase):
    def setUp(self):
        self.record = DummyRecord(seq="A"*300)
        self.cdses = {
            "head": DummyCDS(locus_tag="head", start=0, end=90),
            "tail": DummyCDS(locus_tag="tail", start=200, end=290),
        }
        for cds in self.cdses.values():
            self.record.add_cds_feature(cds)

        self.domains = {
            "head": [DummyHMMResult("PKS_KS", start=0, end=10), DummyHMMResult("PKS_AT", start=20, end=30)],
            "tail": [DummyHMMResult("PKS_ER", start=0, end=10), DummyHMMResult("PP-binding", start=10, end=20)],
        }
        self.ks_subtypes = {'head': [DummyHMMResult("Iterative-KS", start=0, end=10)]}

    def test_merge(self, _patched_fasta, _patched_motifs):
        cdses = self.record.get_cds_features()
        with patch.object(Record, "get_cds_features_within_regions", return_value=cdses):
            with patch.object(domain_identification, "find_domains", return_value=self.domains):
                with patch.object(domain_identification, "find_subtypes", return_value=self.ks_subtypes):
                    results = domain_identification.generate_domains(self.record)
        assert results
        head = results.cds_results[self.cdses["head"]]
        tail = results.cds_results[self.cdses["tail"]]
        assert len(head.modules) == 1
        assert len(tail.modules) == 0  # should have merged into the head
        module = head.modules[0]
        assert module.is_complete()  # only because it was merged

    def test_interrupted_merge(self, _patched_fasta, _patched_motifs):
        # insert the interrupt CDS between the two "halves" of the potential module
        spacer_cds = DummyCDS(locus_tag="mid", start=100, end=130)
        self.record.add_cds_feature(spacer_cds)

        cdses = self.record.get_cds_features()
        with patch.object(Record, "get_cds_features_within_regions", return_value=cdses):
            with patch.object(domain_identification, "find_domains", return_value=self.domains):
                with patch.object(domain_identification, "find_subtypes", return_value=self.ks_subtypes):
                    results = domain_identification.generate_domains(self.record)
        assert results
        head = results.cds_results[self.cdses["head"]]
        assert spacer_cds not in results.cds_results
        tail = results.cds_results[self.cdses["tail"]]
        assert len(head.modules) == 1
        assert len(tail.modules) == 1  # should *not* have merged into head
        assert not head.modules[0].is_complete()
        assert not tail.modules[0].is_complete()
