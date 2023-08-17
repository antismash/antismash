# License: GNU Affero General Public License v3 or later
# A copy of GNU AGPL v3 should have been included in this software package in LICENSE.txt.

# for test files, silence irrelevant and noisy pylint warnings
# pylint: disable=use-implicit-booleaness-not-comparison,protected-access,missing-docstring

from copy import deepcopy
import glob
import os
import unittest
from unittest.mock import patch

from helperlibs.wrappers.io import TemporaryDirectory

from antismash import main
from antismash.common import path, secmet
from antismash.common.module_results import ModuleResults
from antismash.common.test import helpers
from antismash.common.subprocessing.diamond import run_diamond_version
from antismash.config import build_config, get_config, update_config, destroy_config
from antismash.detection import hmm_detection
from antismash.modules import clusterblast
from antismash.modules.clusterblast import core, known
from antismash.outputs import html


def known_dir(filename, *_args):
    return path.get_full_path(__file__, "data", "known", filename)


class Base(unittest.TestCase):
    def setUp(self):
        options = build_config(self.get_args(), isolated=True, modules=[clusterblast, html, hmm_detection])
        _major, _minor, _patch = map(int, run_diamond_version().split("."))
        self.diamond_ver_major = _major
        self.diamond_ver_minor = _minor
        self.diamond_ver_patch = _patch
        self.diamond_version = (_major, _minor, _patch)
        self.old_config = get_config().__dict__
        update_config({"genefinding_gff3": ""})
        self.options = update_config(options)

        with patch.object(known, "_get_datafile_path", side_effect=known_dir):
            assert clusterblast.check_prereqs(self.options) == []
        assert clusterblast.check_options(self.options) == []
        assert clusterblast.is_enabled(self.options)

    def get_args(self):
        """ override with the args you'll need to use in setUp(),
            format is as on the commandline, e.g. ["--tta", "--minimal"]
        """
        raise NotImplementedError("get_args not overridden")

    def tearDown(self):
        destroy_config()
        update_config(self.old_config)

    def get_results(self, results):
        """ override with a function that fetches specific results instance
            for further testing as a tuple of
            (specific_results, global_results)
        """
        raise NotImplementedError("get_results not overridden")

    def check_svgs(self, results, expected, svg_dir):
        # make sure no svgs created yet
        assert not glob.glob(os.path.join(svg_dir, "*.svg"))
        results.write_svg_files(svg_dir)
        # check there's an svg for each result (up to the limit) + 1 combined
        num_svgs = len(glob.glob(os.path.join(svg_dir, "*.svg")))
        assert num_svgs == min(self.options.cb_nclusters, expected) + 1

    def run_antismash(self, filename, expected):
        annotated_records = []

        def callback(tempdir):
            output_gbk = os.path.join(tempdir, os.path.basename(filename))
            annotated_records.extend(secmet.Record.from_genbank(output_gbk))

        with TemporaryDirectory() as output_dir:
            update_config({"output_dir": output_dir})
            with patch.object(main, "get_all_modules", return_value=[hmm_detection, clusterblast, html]):
                with patch.object(main, "_get_all_enabled_modules", return_value=[hmm_detection, clusterblast, html]):
                    run_and_regen = helpers.run_and_regenerate_results_for_module
                    results = run_and_regen(filename, clusterblast, self.options, callback=callback)
            assert annotated_records
            update_config({"output_dir": ""})
            results, global_results = self.get_results(results)
            assert len(results.region_results) == 1
            cluster = results.region_results[0]
            assert len(cluster.ranking) == expected  # will change if database does
            self.check_svgs(global_results, expected, output_dir)
        return annotated_records, results

    def check_nisin(self, expected):
        return self.run_antismash(helpers.get_path_to_nisin_genbank(), expected)

    def check_balhymicin(self, expected):
        return self.run_antismash(helpers.get_path_to_balhymicin_genbank(), expected)


@patch.object(known, "_get_datafile_path", side_effect=known_dir)
class GeneralIntegrationTest(Base):
    def get_args(self):
        return ["--cb-general", "--minimal", "--enable-html",
                "--data", path.get_full_path(__file__, "data")]

    def get_results(self, results):
        assert isinstance(results, ModuleResults)
        assert results.general
        assert results.general.search_type == "clusterblast"
        assert results.knowncluster is None
        assert results.subcluster is None
        return results.general, results

    def test_nisin(self, _patched_known):
        expected_hits = 3 if self.diamond_version < (2, 0, 15) else 2
        records, results = self.check_nisin(expected_hits)
        assert len(records) == 1
        record = records[0]
        ranking = results.region_results[0].ranking

        # check the JSON round trip
        raw = results.to_json()
        # check some values
        for _, score in raw["results"][0]["ranking"]:
            assert score["similarity"] > 0
        # and check it parses back correctly
        rebuilt = results.from_json(deepcopy(raw), record)
        assert results.region_results[0].ranking[0][1].hits == 11
        assert isinstance(rebuilt, type(results))
        new_json = rebuilt.to_json()
        assert new_json["results"][0]["ranking"]
        for _, score in new_json["results"][0]["ranking"]:
            assert score["similarity"] > 0
        assert new_json == raw

        assert len(ranking) == expected_hits
        match, score = ranking[0]
        assert match.accession == "NC_017486"
        assert match.cluster_label == "c1234-56789"
        assert score.score == 30
        assert score.hits == 11
        assert score.core_gene_hits == 2
        assert score.blast_score == 8401.0 if self.diamond_ver_major < 2 else 8402.0
        assert score.synteny_score == 14
        assert score.core_bonus == 3
        assert len(score.scored_pairings) == score.hits
        query, subject = score.scored_pairings[0]
        assert query.id == "nisA"
        assert subject.name == "CVCAS_RS03115"

        match, score = ranking[1]
        assert match.accession == "ALTERED_NC_017486"
        assert match.cluster_label == "c1234-56789"
        assert score.score == 24
        assert score.hits == 10
        assert score.core_gene_hits == 1
        assert score.blast_score == 7579.0 if self.diamond_ver_major < 2 else 7577.0
        assert score.synteny_score == 10
        assert score.core_bonus == 3
        assert len(score.scored_pairings) == score.hits
        query, subject = score.scored_pairings[0]
        assert query.id == "nisA"
        assert subject.name == "CVCAS_RS03115"


@patch.object(known, "_get_datafile_path", side_effect=known_dir)
class KnownIntegrationTest(Base):
    def get_args(self):
        return ["--cb-knowncluster", "--minimal", "--enable-html"]

    def get_results(self, results):
        assert isinstance(results, ModuleResults)
        assert results.general is None
        assert results.knowncluster
        assert results.knowncluster.search_type == "knownclusterblast"
        assert results.subcluster is None
        return results.knowncluster, results

    def test_nisin(self, _patched_dir):
        # blast scores not checked due to diamond versions having different results
        _, results = self.check_nisin(2)
        ranking = results.region_results[0].ranking
        assert len(ranking) == 2
        match, score = ranking[0]
        assert match.accession == "BGC0000536"
        assert match.cluster_label == "c1"
        assert score.score == 30
        assert score.hits == 11
        assert score.core_gene_hits == 2
        assert score.synteny_score == 14
        assert score.core_bonus == 3
        assert len(score.scored_pairings) == score.hits
        query, subject = score.scored_pairings[0]
        assert query.id == "nisA"
        assert subject.name == "BAG71479.1"

        match, score = ranking[1]
        assert match.accession == "BGC0000549"
        assert match.cluster_label == "c1"
        assert score.score == 23
        assert score.hits == 9
        assert score.core_gene_hits == 2
        assert score.synteny_score == 9
        assert score.core_bonus == 3
        assert len(score.scored_pairings) == 10  # some having shared hits
        query, subject = score.scored_pairings[0]
        assert query.id == "nisA"
        assert subject.name == "AEX55166.1"

    def test_balhymicin(self, _patched_dir):
        # blast scores and derivative values not checked due to diamond versions having different results
        _, results = self.check_balhymicin(2)
        ranking = results.region_results[0].ranking
        assert len(ranking) == 2
        match, score = ranking[0]
        assert match.accession == "BGC0000455"
        assert match.cluster_label == "c1"
        assert score.hits == 31
        assert score.core_gene_hits == 5
        assert score.core_bonus == 3
        assert len(score.scored_pairings) >= 31  # some having shared hits
        query, subject = score.scored_pairings[0]
        assert query.id == "bbr"
        assert subject.name == "AEI58862.1"

        match, score = ranking[1]
        assert match.accession == "BGC0001178"
        assert match.cluster_label == "c1"
        assert score.hits == 29
        assert score.core_gene_hits == 5
        assert score.core_bonus == 3
        assert len(score.scored_pairings) >= score.hits
        query, subject = score.scored_pairings[0]
        assert query.id == "vanS"
        assert subject.name == "AGS77301.1"

    def test_fusariam_scirpi(self, _patched_dir):
        # this is a special case where it's a single CDS cluster that matches
        # against other single CDS clusters, before this test they were all
        # discarded
        genbank = path.get_full_path(__file__, "data", "Z18755.3.gbk")
        _, results = self.run_antismash(genbank, 2)
        assert list(results.mibig_entries) == [1]  # only one region in record
        assert list(results.mibig_entries[1]) == ["esyn1"]  # and only one CDS
        # 2 hits against single-CDS MiBIG clusters, only those are reported
        # as multi-CDS reference clusters need multiple query CDSs to hit
        for ref_cluster, _ in results.region_results[0].ranking:
            assert len(ref_cluster.proteins) == 1


class SubIntegrationTest(Base):
    def get_args(self):
        return ["--cb-subcluster", "--minimal", "--enable-html"]

    def get_results(self, results):
        assert isinstance(results, ModuleResults)
        assert results.general is None
        assert results.knowncluster is None
        assert results.subcluster
        assert results.subcluster.search_type == "subclusterblast"
        return results.subcluster, results

    def test_nisin(self):
        self.check_nisin(0)

    def test_balhymicin(self):
        self.check_balhymicin(21)


class TestDatabaseValidity(unittest.TestCase):
    def setUp(self):
        options = build_config([], isolated=True, modules=[clusterblast])
        self.old_config = get_config().__dict__
        self.options = update_config(options)

    def tearDown(self):
        destroy_config()
        update_config(self.old_config)

    def _check_proteins_match_clusters(self, searchtype):
        clusters = core.load_reference_clusters(searchtype)
        proteins = core.load_reference_proteins(searchtype)
        for cluster in clusters.values():
            for protein in cluster.tags:
                assert protein in proteins, f"missing: {protein}"

    def test_general(self):
        self._check_proteins_match_clusters("clusterblast")

    def test_known(self):
        self._check_proteins_match_clusters("knownclusterblast")

    def test_sub(self):
        self._check_proteins_match_clusters("subclusterblast")
