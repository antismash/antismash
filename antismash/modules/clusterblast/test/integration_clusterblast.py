# License: GNU Affero General Public License v3 or later
# A copy of GNU AGPL v3 should have been included in this software package in LICENSE.txt.

# for test files, silence irrelevant and noisy pylint warnings
# pylint: disable=no-self-use,protected-access,missing-docstring

import glob
import os
import unittest
from unittest.mock import patch

from helperlibs.wrappers.io import TemporaryDirectory

from antismash.main import get_all_modules
from antismash.common import path, subprocessing
from antismash.common.module_results import ModuleResults
import antismash.common.test.helpers as helpers
from antismash.config import build_config, get_config, update_config, destroy_config
from antismash.modules import clusterblast
from antismash.modules.clusterblast import core, known

MOCKED_DATA = path.get_full_path(__file__, "data")


class Base(unittest.TestCase):
    def setUp(self):
        options = build_config(self.get_args(), isolated=True, modules=get_all_modules())
        self.old_config = get_config().__dict__
        self.options = update_config(options)

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
        with TemporaryDirectory() as output_dir:
            update_config({"output_dir": output_dir})
            results = helpers.run_and_regenerate_results_for_module(filename, clusterblast, self.options)
            update_config({"output_dir": ""})
            results, global_results = self.get_results(results)
            assert len(results.region_results) == 1
            cluster = results.region_results[0]
            assert len(cluster.ranking) == expected  # will change if database does
            self.check_svgs(global_results, expected, output_dir)
        return results

    def check_nisin(self, expected):
        return self.run_antismash(helpers.get_path_to_nisin_genbank(), expected)

    def check_balhymicin(self, expected):
        return self.run_antismash(helpers.get_path_to_balhymicin_genbank(), expected)


class GeneralIntegrationTest(Base):
    def get_args(self):
        return ["--cb-general", "--minimal", "--data", path.get_full_path(__file__, "data")]

    def get_results(self, results):
        assert isinstance(results, ModuleResults)
        assert results.general
        assert results.general.search_type == "clusterblast"
        assert results.knowncluster is None
        assert results.subcluster is None
        return results.general, results

    def test_nisin(self):
        results = self.check_nisin(2)
        ranking = results.region_results[0].ranking
        assert len(ranking) == 2
        match, score = ranking[0]
        assert match.accession == "NC_017486"
        assert match.cluster_label == "c2"
        assert score.score == 30
        assert score.hits == 11
        assert score.core_gene_hits == 2
        assert score.blast_score == 8401.0
        assert score.synteny_score == 14
        assert score.core_bonus == 3
        assert len(score.scored_pairings) == score.hits
        query, subject = score.scored_pairings[0]
        assert query.id == "nisA"
        assert subject.name == "CVCAS_RS03115"

        match, score = ranking[1]
        assert match.accession == "ALTERED_NC_017486"
        assert match.cluster_label == "c2"
        assert score.score == 24
        assert score.hits == 10
        assert score.core_gene_hits == 1
        assert score.blast_score == 7579.0
        assert score.synteny_score == 10
        assert score.core_bonus == 3
        assert len(score.scored_pairings) == score.hits
        query, subject = score.scored_pairings[0]
        assert query.id == "nisA"
        assert subject.name == "CVCAS_RS03115"


@patch.object(known, "_SHIPPED_DATA_DIR", MOCKED_DATA)
class KnownIntegrationTest(Base):
    def get_args(self):
        return ["--cb-knowncluster", "--minimal"]

    def get_results(self, results):
        assert isinstance(results, ModuleResults)
        assert results.general is None
        assert results.knowncluster
        assert results.knowncluster.search_type == "knownclusterblast"
        assert results.subcluster is None
        return results.knowncluster, results

    def test_nisin(self):
        # blast scores not checked due to diamond versions having different results
        results = self.check_nisin(2)
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

    def test_balhymicin(self):
        # blast scores and derivative values not checked due to diamond versions having different results
        results = self.check_balhymicin(2)
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

    def test_fusariam_scirpi(self):
        # this is a special case where it's a single CDS cluster that matches
        # against other single CDS clusters, before this test they were all
        # discarded
        genbank = path.get_full_path(__file__, "data", "Z18755.3.gbk")
        results = self.run_antismash(genbank, 2)
        assert list(results.mibig_entries) == [1]  # only one region in record
        assert list(results.mibig_entries[1]) == ["esyn1"]  # and only one CDS
        # 2 hits against single-CDS MiBIG clusters, only those are reported
        # as multi-CDS reference clusters need multiple query CDSs to hit
        for ref_cluster, _ in results.region_results[0].ranking:
            assert len(ref_cluster.proteins) == 1


class SubIntegrationTest(Base):
    def get_args(self):
        return ["--cb-subcluster", "--minimal"]

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


class TestDiamondDatabaseChecks(unittest.TestCase):
    def setUp(self):
        self.format0_file = path.get_full_path(__file__, "data", "format0.dmnd")
        self.format1_file = path.get_full_path(__file__, "data", "format1.dmnd")
        self.empty = path.get_full_path(__file__, "data", "empty.dmnd")

        options = build_config([], isolated=True, modules=get_all_modules())
        self.old_config = get_config().__dict__
        self.options = update_config(options)

    def tearDown(self):
        destroy_config()
        update_config(self.old_config)

    def test_check_diamond_db_compatible(self):
        with TemporaryDirectory(change=True):
            dummy_fasta = "dummy.fa"
            dummy_db = "dummy.dmnd"
            with open(dummy_fasta, "w") as handle:
                handle.write(">test\nM\n")
            subprocessing.run_diamond_makedb(dummy_db, dummy_fasta)
            compatible_format = core._extract_db_format(dummy_db)
            assert core.check_diamond_db_compatible(dummy_db)

        broken_file = self.format0_file if compatible_format > 0 else self.format1_file

        assert not core.check_diamond_db_compatible(broken_file)
        assert not core.check_diamond_db_compatible(self.empty)


class TestDatabaseValidity(unittest.TestCase):
    def setUp(self):
        options = build_config([], isolated=True, modules=get_all_modules())
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
                assert protein in proteins, "missing: %s" % protein

    def test_general(self):
        self._check_proteins_match_clusters("clusterblast")

    def test_known(self):
        self._check_proteins_match_clusters("knownclusterblast")

    def test_sub(self):
        self._check_proteins_match_clusters("subclusterblast")
