# License: GNU Affero General Public License v3 or later
# A copy of GNU AGPL v3 should have been included in this software package in LICENSE.txt.

# for test files, silence irrelevant and noisy pylint warnings
# pylint: disable=no-self-use,protected-access,missing-docstring

import glob
import os
import unittest

from helperlibs.wrappers.io import TemporaryDirectory

from antismash.main import get_all_modules
from antismash.common import path
from antismash.common.module_results import ModuleResults
import antismash.common.test.helpers as helpers
from antismash.config import build_config, get_config, update_config, destroy_config
from antismash.modules import clusterblast


class Base(unittest.TestCase):
    def setUp(self):
        options = build_config(self.get_args(), isolated=True, modules=get_all_modules())
        self.old_config = get_config().__dict__
        self.options = update_config(options)

        assert clusterblast.check_prereqs() == []
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
            if expected > clusterblast.get_result_limit():
                assert cluster.total_hits == expected
                expected = clusterblast.get_result_limit()
            assert len(cluster.ranking) == expected  # will change if database does
            self.check_svgs(global_results, expected, output_dir)
        return results

    def check_nisin(self, expected):
        return self.run_antismash(helpers.get_path_to_nisin_genbank(), expected)

    def check_balhymicin(self, expected):
        return self.run_antismash(helpers.get_path_to_balhymicin_genbank(), expected)


class GeneralIntegrationTest(Base):
    # TODO: test with a small sequence instead (grab a CDS that hit and take it's translation)
    def get_args(self):
        return ["--cb-general", "--minimal"]

    def get_results(self, results):
        assert isinstance(results, ModuleResults)
        assert results.general
        assert results.general.search_type == "clusterblast"
        assert results.knowncluster is None
        assert results.subcluster is None
        return results.general, results

    def test_nisin(self):
        self.check_nisin(2452)


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
        self.check_nisin(7)

    def test_balhymicin(self):
        self.check_balhymicin(102)

    def test_fusariam_scirpi(self):
        # this is a special case where it's a single CDS cluster that matches
        # against other single CDS clusters, before this test they were all
        # discarded
        genbank = path.get_full_path(__file__, "data", "Z18755.3.gbk")
        results = self.run_antismash(genbank, 2)
        assert list(results.mibig_entries) == [1]  # only one region in record
        assert list(results.mibig_entries[1]) == ["CAA79245.2"]  # and only one CDS
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
