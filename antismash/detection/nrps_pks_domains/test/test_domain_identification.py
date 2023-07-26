# License: GNU Affero General Public License v3 or later
# A copy of GNU AGPL v3 should have been included in this software package in LICENSE.txt.

# for test files, silence irrelevant and noisy pylint warnings
# pylint: disable=use-implicit-booleaness-not-comparison,protected-access,missing-docstring

import os
import unittest
from unittest.mock import patch

from antismash.common import fasta, path, subprocessing, utils
from antismash.common.secmet import Record
from antismash.common.test.helpers import (
    DummyCDS,
    DummyHMMResult,
    DummyRecord,
    DummyRegion,
    get_simple_options,
)
from antismash.common.secmet.test.helpers import (
    DummySubRegion,
)
from antismash.config import destroy_config, update_config
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
                    assert val == 1, f"{ks_name} {vars(counter)}"
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
@patch.object(domain_identification, "get_database_path", return_value="")
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

    def test_merge(self, _patched_fasta, _patched_motifs, _patched_dbs):
        cdses = self.record.get_cds_features()
        with patch.object(Record, "get_cds_features_within_regions", return_value=cdses):
            with patch.object(domain_identification, "find_domains", return_value=self.domains):
                with patch.object(domain_identification, "find_subtypes", side_effect=[self.ks_subtypes, []]):
                    results = domain_identification.generate_domains(self.record)
        assert results
        head = results.cds_results[self.cdses["head"]]
        tail = results.cds_results[self.cdses["tail"]]
        assert len(head.modules) == 1
        assert len(tail.modules) == 0  # should have merged into the head
        module = head.modules[0]
        assert module.is_complete()  # only because it was merged

    def test_interrupted_merge(self, _patched_fasta, _patched_motifs, _patched_dbs):
        # insert the interrupt CDS between the two "halves" of the potential module
        spacer_cds = DummyCDS(locus_tag="mid", start=100, end=130)
        self.record.add_cds_feature(spacer_cds)

        cdses = self.record.get_cds_features()
        with patch.object(Record, "get_cds_features_within_regions", return_value=cdses):
            with patch.object(domain_identification, "find_domains", return_value=self.domains):
                with patch.object(domain_identification, "find_subtypes", side_effect=[self.ks_subtypes, []]):
                    results = domain_identification.generate_domains(self.record)
        assert results
        head = results.cds_results[self.cdses["head"]]
        assert spacer_cds not in results.cds_results
        tail = results.cds_results[self.cdses["tail"]]
        assert len(head.modules) == 1
        assert len(tail.modules) == 1  # should *not* have merged into head
        assert not head.modules[0].is_complete()
        assert not tail.modules[0].is_complete()

    def test_merging_across_regions(self, _patched_fasta, _patched_motifs, _patched_dbs):
        # make sure the two halves are in different regions
        for cds in self.cdses.values():
            subregion = DummySubRegion(start=cds.location.start, end=cds.location.end)
            self.record.add_subregion(subregion)
            region = DummyRegion(candidate_clusters=[], subregions=[subregion])
            self.record.add_region(region)
            assert cds.region == region

        with patch.object(domain_identification, "find_domains", return_value=self.domains):
            with patch.object(domain_identification, "find_subtypes", side_effect=[self.ks_subtypes, []]):
                results = domain_identification.generate_domains(self.record)
        # the two should not have been merged
        for cds in self.cdses.values():
            modules = results.cds_results[cds].modules
            assert len(modules) == 1
            assert not modules[0].is_complete()


class TestDatabases(unittest.TestCase):
    def tearDown(self):
        domain_identification.DATABASE_PATHS.clear()
        destroy_config()

    def test_caching(self):
        root = "/some/dummy/root"
        update_config(get_simple_options(None, ['--databases', root]))
        assert not domain_identification.DATABASE_PATHS
        with patch.object(path, "find_latest_database_version",
                          return_value="1.0") as patched:
            result = domain_identification.get_database_path("some_dir", "some.file")
            assert patched.called
        expected_subdir = os.path.join(root, "nrps_pks", "some_dir", "1.0")
        expected_full_path = os.path.join(expected_subdir, "some.file")
        assert result == expected_full_path
        assert domain_identification.DATABASE_PATHS == {"some_dir": expected_subdir}

        # check that the find isn't used while we build paths for different files
        with patch.object(path, "find_latest_database_version",
                          return_value="bad") as patched:
            for filename in ["some.file", "other"]:
                result = domain_identification.get_database_path("some_dir", filename)
                assert not patched.called
                assert result == os.path.join(expected_subdir, filename)


class TestTransatorPrediction(unittest.TestCase):
    def setUp(self):
        self.record = DummyRecord(seq="A"*20)
        self.cds = DummyCDS(start=1, end=20)
        self.record.add_cds_feature(self.cds)

    def build_hit_from_types(self, types: list[str]) -> DummyHMMResult:
        types = list(types)
        start = 1
        end = len(types) * 2
        current = DummyHMMResult(types.pop(), start=start, end=end)
        for outer in reversed(types):
            end -= 1
            start += 1
            current = DummyHMMResult(outer, start=start, end=end, internal_hits=[current])
        return current

    def annotate(self, types: list[str]) -> None:
        hit = self.build_hit_from_types(types)
        assert hit.detailed_names == types
        result = domain_identification.CDSResult([hit], [], [])
        result.annotate_domains(self.record, self.cds)

    def test_irrelevant(self):
        types = ["PKS_KS", "some subtype", "innermost"]
        self.annotate(types)
        assert self.cds.nrps_pks.domains[0].subtypes == types[1:]
        assert not self.cds.nrps_pks.domains[0].get_predictions()

    def test_with_prediction(self):
        types = ["PKS_KS", "Trans-AT-KS", "innermost"]
        self.annotate(types)
        assert self.cds.nrps_pks.domains[0].subtypes == types[1:]
        assert self.cds.nrps_pks.domains[0].get_predictions()["transATor"] == types[-1]

    def test_without_prediction(self):
        types = ["PKS_KS", "Trans-AT-KS"]
        self.annotate(types)
        assert self.cds.nrps_pks.domains[0].subtypes == types[1:]
        assert self.cds.nrps_pks.domains[0].get_predictions()["transATor"] == "(unknown)"
